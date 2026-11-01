//! Native serial caller, tested through the internal probe before public dispatch.
#[cfg(feature = "parallel")]
use super::reader_batch::{
    attach, prepare_lone_unit, prepare_reader_batch, prepare_unit, ReaderBatchShape, UnitPrep,
};
use super::{
    archive,
    args::{BlastxArgs, SemanticOptionsError},
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    // ```c++
    // bool CStreamLineReader::AtEOF(void) const
    // {
    //     return !m_UngetLine &&
    //         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
    // }
    // ```
    input::{
        parse_fasta_stream_keeping, stream_is_empty, FastaParseError, FastaStream, InputWarning,
    },
    query_setup::BATCH_SIZE,
    report,
    runtime::search_internal,
};
use super::{
    database::{BlastxDatabase, ReportedSubjects},
    input::FastaRecord,
    report::ShownSubjects,
    runtime::BatchResults,
    subjects::{Subjects, UnreadableSubjects},
};
use crate::algorithm::common::database::{
    prefetches_subjects, HeaderReadError, PreliminaryReadError, TaxonomyRefusal,
};
use crate::report::blast_xml2::{FsSideFiles, NoSideFiles, RunSideFiles, Xml2Kind};
use anyhow::{bail, Result};
use std::{
    fmt,
    fs::File,
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    // ```c++
    // bool CStreamLineReader::AtEOF(void) const
    // {
    //     return !m_UngetLine &&
    //         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
    // }
    // ```
    io::{self, BufWriter, Write},
    path::Path,
};

// The NCBI-style exit error is shared with the other programs' command lines.
use crate::cli::inaccessible;
pub use crate::cli::{exit_on_native_error, NativeError};
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
// ```c++
//         } else {                                                            \
//             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
//             exit_code = BLAST_ENGINE_ERROR;                                 \
// ```
fn engine(message: impl fmt::Display) -> anyhow::Error {
    NativeError {
        exit: 3,
        message: format!("BLAST engine error: {message}\n"),
    }
    .into()
}

/// The error of the database's sequences (`BlastxDatabase::x_check_stores`, or a read of
/// `x_db_stream.rs`): NCBI's engine error of a reader error, or LOSAT's refusal as it is.
fn database_error(kept: anyhow::Error) -> anyhow::Error {
    match kept.downcast::<PreliminaryReadError>() {
        Ok(PreliminaryReadError(text)) => engine(text),
        Err(refusal) => refusal,
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:181-184
// ```c++
//     catch (const CObjReaderParseException& e) {                             \
//         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
fn input_error(error: anyhow::Error) -> anyhow::Error {
    if let Some(e) = error.downcast_ref::<FastaParseError>() {
        NativeError {
            exit: 1,
            message: format!("BLAST query error: {e}\n"),
        }
        .into()
    } else {
        error
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
// ```c++
//     }                                                                       \
//     catch (const std::ios::failure&) {                                      \
//         LOG_POST(Error << "BLAST failed to write output");                  \
//         exit_code = BLAST_OUTPUT_ERROR;                                     \
//     }                                                                       \
// ```
fn output_error(_: io::Error) -> anyhow::Error {
    NativeError {
        exit: 6,
        message: "BLAST failed to write output\n".into(),
    }
    .into()
}
/// The error of writing one batch's archive: the encoder's refusal of a byte NCBI writes as `#`
/// (`SerialError`, after the query titles were checked before the search; reachable only when
/// the query file changes between that read and the search) is LOSAT's refusal, not an error of
/// NCBI's engine, so it gets no `BLAST engine error:` and no error archive (`archive::error_lines`
/// takes `NativeError`s only), as `XmlRefusal` gets none from `report_error`.
fn archive_error(error: anyhow::Error) -> anyhow::Error {
    if error.is::<crate::report::serial::encoder::SerialError>() {
        error
    } else {
        report_error(error)
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
// ```c++
//         } else {                                                            \
//             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
//             exit_code = BLAST_ENGINE_ERROR;                                 \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
// ```c++
//     }                                                                       \
//     catch (const std::ios::failure&) {                                      \
//         LOG_POST(Error << "BLAST failed to write output");                  \
//         exit_code = BLAST_OUTPUT_ERROR;                                     \
//     }                                                                       \
// ```
fn report_error(error: anyhow::Error) -> anyhow::Error {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
    // ```c++
    // CBlastTabularInfo::~CBlastTabularInfo()
    // {
    //     m_Ostream.flush();
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
    // ```c++
    // {
    //     m_Outfile.exceptions(NcbiBadbit);
    // ```
    if error.is::<super::blast_xml::XmlRefusal>() {
        // LOSAT's refusal of a title NCBI's XML cannot write, not an error of NCBI's engine.
        error
    } else if error.downcast_ref::<NativeError>().is_some_and(|e| {
        (e.exit == 1 && e.message.starts_with("Command line argument error"))
            || (e.exit == 6 && e.message == "BLAST failed to write output\n")
    }) {
        // An outfmt 13 or 14 side file NCBI cannot open (`CArgException`, exit status 1), or a
        // write to one that fails (`framing::side_file_write_error`).
        error
    } else if error.downcast_ref::<report::TabularFlushError>().is_some() {
        NativeError{exit:-6,message:"terminate called after throwing an instance of 'std::__ios_failure'\n  what():  basic_ios::clear: iostream error\n".into()}.into()
    } else if let Some(e) = error.downcast_ref::<io::Error>() {
        output_error(io::Error::new(e.kind(), e.to_string()))
    } else {
        engine(error)
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:62-70
// ```c++
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
//     m_BlastDbArgs->SetDatabaseMaskingSupport(true);
//     m_BlastDbArgs->SetIPGFilteringSupport(true);
//     arg.Reset(m_BlastDbArgs);
//     m_Args.push_back(arg);
//
//     m_StdCmdLineArgs.Reset(new CStdCmdLineArgs);
//     arg.Reset(m_StdCmdLineArgs);
//     m_Args.push_back(arg);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2537-2537
// ```c++
//             subj_input_stream = &args[kArgSubject].AsInputFile();
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2554-2557
// ```c++
//         m_Scope = ReadSequencesToBlast(*subj_input_stream, IsProtein(),
//                                        subj_range, parse_deflines,
//                                        use_lcase_masks, subjects, m_IsMapper);
//         m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects));
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3468-3468
// ```c++
//             m_InputStream = &args[kArgQuery].AsInputFile();
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3479-3479
// ```c++
//         m_OutputStream = &args[kArgOutput].AsOutputFile();
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:198-200
// ```c++
//         CRef<CScope> scope;
//         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
//                          db_adapter, scope);
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:213-217
// ```c++
//         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
//            	ERR_POST(Warning << "Query is Empty!");
//            	return BLAST_EXIT_SUCCESS;
//         }
//         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
// ```c++
//         formatter.PrintProlog();
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:260-297
// ```c++
// 	    BLAST_PROF_START( APP.LOOP.PRE );
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//
//             SaveSearchStrategy(args, m_CmdLineArgs, queries, m_OptsHndl);
//
//             CRef<CSearchResultSet> results;
// 	    BLAST_PROF_STOP( APP.LOOP.PRE );
//
//             if (m_CmdLineArgs->ExecuteRemotely()) {
//                 CRef<CRemoteBlast> rmt_blast =
//                     InitializeRemoteBlast(queries, db_args, m_OptsHndl,
//                           m_CmdLineArgs->ProduceDebugRemoteOutput(),
//                           m_CmdLineArgs->GetClientId());
//                 results = rmt_blast->GetResultSet();
//             } else {
// 	        BLAST_PROF_START( APP.LOOP.BLAST );
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
// 	        BLAST_PROF_STOP( APP.LOOP.BLAST );
//             }
// 	    BLAST_PROF_START( APP.LOOP.FMT );
//             if (fmt_args->ArchiveFormatRequested(args)) {
//                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
//                 m_Bah.ResetMessages();
//             } else {
//                 BlastFormatter_PreFetchSequenceData(*results, scope,
//                 		                            fmt_args->GetFormattedOutputChoice());
//             	ITERATE(CSearchResultSet, result, *results) {
//                	    formatter.PrintOneResultSet(**result, query_batch);
//             	}
//             }
// 	    BLAST_PROF_STOP( APP.LOOP.FMT );
// 	    batch_num++;
//         }
//         BLAST_PROF_START( APP.POST );
//         formatter.PrintEpilog(opt);
// ```
fn run_inner(
    args: &BlastxArgs,
    output: &mut Option<BufWriter<Box<dyn Write>>>,
    archive_slot: &mut Option<crate::report::archive::ArchiveHandle>,
) -> Result<()> {
    // Diagnostic (LOSAT_TIMING): the clock of the run (x_timing.rs). No NCBI counterpart.
    let x_t_run = super::x_timing::lap();
    use super::x_timing::Run as XRun;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3631-3639
    // ```c++
    //     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
    //         (*arg)->ExtractAlgorithmOptions(args, opts);
    //     }
    //
    //     m_IsUngapped = !opts.GetGappedMode();
    //     try { retval->Validate(); }
    //     catch (const CBlastException& e) {
    //         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
    // ```c++
    //     if (GetExportSearchStrategyStream(args) ||
    //            m_FormattingArgs->ArchiveFormatRequested(args)) {
    //         locality = CBlastOptions::eBoth;
    //     }
    // ```
    // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
    // option handlers run, with NCBI's errors; a format that NCBI runs and LOSAT does not
    // write is rejected there too. The custom tabular fields are read as NCBI reads them
    // (`ResolvedOptions::tabular`).
    use crate::blastinput::app::{self, Program};
    let choice = app::parse_formatting_string(&args.outfmt)?;
    let format = app::report_format(&choice, Program::Blastx, args.out.as_deref())?;
    // The BLAST archive (outfmt 11): its request, and the refusals of what LOSAT cannot write,
    // before anything is read or written.
    // The error archives follow the `-outfmt` text (`writes_error_archive`).
    let archive = archive::handle_for(args, format == Some(app::ReportFormat::ArchiveFormat))?
        .map(|handle| handle.with_outfmt_text(&args.outfmt));
    archive_slot.clone_from(&archive);
    // Preserve explicit product refusals; defer supported semantic failures
    // until subject and standard query/output streams have been extracted.
    let resolution = args.resolve();
    if resolution
        .as_ref()
        .is_err_and(|e| e.downcast_ref::<SemanticOptionsError>().is_none())
    {
        return resolution.map(|_| ());
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3157-3165
    // ```c++
    //
    //     arg_desc.AddDefaultKey(kArgNumThreads, "int_value",
    //                            "Number of threads (CPUs) to use in the BLAST search",
    //                            CArgDescriptions::eInteger,
    //                            NStr::IntToString(kDfltValue));
    //     arg_desc.SetConstraint(kArgNumThreads,
    //                            new CArgAllowValuesGreaterThanOrEqual(kMinValue));
    //     arg_desc.SetDependency(kArgNumThreads,
    //                            CArgDescriptions::eExcludes,
    // ```
    // Existing LOSAT search-scoped pool validates the explicit thread contract
    // before input/output side effects; serial WASI remains explicitly serial.
    crate::utils::threading::validate_threads(args.num_threads as usize)?;
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:727-735
    // ```c++
    //             fstrm->open(AsString().c_str(),IOS_BASE::in | mode);
    //             if ( !fstrm->is_open() ) {
    //                 delete fstrm;
    //                 fstrm = NULL;
    //             } else {
    //                 m_DeleteFlag = true;
    //             }
    //         }
    //         m_Ios = fstrm;
    // ```
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    // ```c++
    // bool CStreamLineReader::AtEOF(void) const
    // {
    //     return !m_UngetLine &&
    //         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
    // }
    // ```
    // `-db`: the database opens after the options are processed (below).
    let subject_path = args.subject.as_deref();
    // NCBI opens the `-subject` and `-query` files inside `Run`, after the diagnostics handler
    // that saves the messages is set: an error archive follows their errors.
    let early_archive_error = |error: anyhow::Error| -> anyhow::Error {
        write_run_error_archive(args, archive.as_ref(), &error);
        error
    };
    let subject_file = match (&args.db, subject_path) {
        (Some(_), _) => None,
        (None, Some(path)) => Some(
            File::open(path)
                .map_err(|_| inaccessible("subject", path))
                .map_err(early_archive_error)?,
        ),
        (None, None) => return Err(crate::blastinput::app::missing_subject_error()),
    };
    // With 2+ threads batches are searched on other threads while this function runs. A diagnostic
    // or a panic message written there must not wait for a lock this thread keeps for the whole
    // run, so the handle is used unlocked: the same bytes, the lock taken for each write.
    // No NCBI counterpart: only the stderr lock is taken per write instead of once for the run, so
    // that other threads can write diagnostics (the bytes and their order are unchanged); it does
    // not change any value NCBI computes.
    // With the BLAST archive the diagnostics also go to its message log (`DiagTee`): the
    // messages of each archive are the diagnostics posted since the previous one.
    let mut archive_stderr = io::stderr();
    let mut archive_tee;
    let mut stderr: Box<dyn Write + '_> = match &archive {
        Some(handle) => {
            archive_tee = crate::report::archive::DiagTee::new(&mut archive_stderr, handle.clone());
            Box::new(&mut archive_tee)
        }
        None => Box::new(io::stderr()),
    };
    // The taxid or seqid list (line DI-8): NCBI's handler of the database arguments
    // (`s_GetTaxIDList`, `CSeqDBFileGiList`) comes before the formatting handler and before
    // `CStdCmdLineArgs`, which opens the query and the output: its warning before theirs (probe
    // F04 of line DI-8), and a list error leaves no `-out`, only the error archive that NCBI's
    // `PrintErrorArchive` writes for outfmt 11 (probe E of line DI-8). NCBI's application layer
    // has read its environment and registry files by then; LOSAT reads them for `-db` (the
    // registry's `BLASTDB`), as its other programs with `-db`, and leaves `-subject`'s path as it
    // was.
    let db_settings = match args.database_name() {
        Some(_) => Some(
            crate::blastinput::ncbi_environment::check_ncbi_application_settings("blastx")
                .map_err(anyhow::Error::msg)?,
        ),
        None => None,
    };
    let list = match &db_settings {
        Some(settings) => args
            .db_filter
            .resolve(settings, &mut stderr)
            .map_err(early_archive_error)?,
        None => None,
    };
    // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let x_t_read_fasta = crate::utils::xtiming::clock();
    let subjects = match subject_file {
        Some(subject_file) => parse_fasta_stream_keeping(
            &mut {
                let mut stream = FastaStream::from_file(subject_file);
                stream.parse_deflines = args.parse_deflines;
                stream
            },
            true,
            args.lcase_masking,
            usize::MAX,
            archive.is_some(),
            &mut |_| Ok(()),
            &mut |warning| reader_warning(warning, &mut stderr),
        )
        .map_err(input_error)
        .map_err(early_archive_error)?,
        None => Vec::new(),
    };
    // LOSAT's refusal of a title the archive cannot write, before its first byte.
    if archive.is_some() {
        archive::check_titles(&subjects, "subject")?;
    }
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, x_t_read_fasta);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    // ```c++
    // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
    //     : m_QueryVector(& queries)
    // {
    //     if (queries.Empty()) {
    //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
    //     }
    // ```
    if args.db.is_none() && subjects.is_empty() {
        return Err(early_archive_error(engine("Empty CBlastQueryVector")));
    }
    let query_file = crate::blastinput::input_files::open_input(&args.query, "query", "BLASTX")
        .map_err(early_archive_error)?;
    let destination: Box<dyn Write> = match args.out.as_deref().filter(|p| *p != Path::new("-")) {
        Some(path) => Box::new(File::create(path).map_err(|_| inaccessible("out", path))?),
        None => Box::new(io::stdout().lock()),
    };
    // BLAST XML (outfmt 5) and the XML2 and JSON documents (outfmt 13 to 16): the output can
    // hold the whole document in memory (`blast_xml::HeldOutput`), decided below before the
    // search.
    // Outfmt 13 and 14 write side files next to `-out` and the master file at `-out`: with no
    // file name they stop at NCBI's `xinclude_check` below (`Please provide a file name`).
    let out_file = args.out.as_deref().filter(|p| *p != Path::new("-"));
    let structured_xml = match choice.number {
        5 | 15 | 16 => true,
        13 | 14 => out_file.is_some(),
        _ => false,
    };
    // The hold of the side files of outfmt 13 and 14 follows the master file's.
    let side_files = out_file.and_then(|out| {
        let kind = app::ReportFormat::from_number(choice.number).and_then(Xml2Kind::of)?;
        kind.is_files()
            .then(|| super::blast_xml::HeldSideFiles::new(FsSideFiles::new(out, kind)))
    });
    let (destination, held) = if structured_xml {
        let held = super::blast_xml::HeldOutput::new(destination, side_files.clone());
        (Box::new(held.clone()) as Box<dyn Write>, Some(held))
    } else {
        (destination, None)
    };
    // The run keeps its output in the caller's slot: an error of NCBI's own is followed by the
    // error archive in the same output (`run`).
    let writer = output.insert(BufWriter::new(destination));
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
    // ```c++
    //     if(hitlist_size < 5){
    //    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886
    // ```c++
    //     ParseFormattingString(args, m_OutputFormat, m_CustomOutputFormatSpec,m_CustomDelim);
    //     if((m_OutputFormat == eSAM) && !(m_FormatFlags & eIsSAM) ){
    //     		NCBI_THROW(CInputException, eInvalidInput,
    //     		                        "SAM format is only applicable to blastn" );
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:110-126
    // ```c++
    //     arg.Reset(new CGeneticCodeArgs(CGeneticCodeArgs::eQuery));
    //     m_Args.push_back(arg);
    //
    //     m_FormattingArgs.Reset(new CFormattingArgs);
    //     arg.Reset(m_FormattingArgs);
    //     m_Args.push_back(arg);
    // ...
    //     arg.Reset(new CCompositionBasedStatsArgs);
    //     m_Args.push_back(arg);
    // ```
    // The formatting handler raises NCBI's errors for the formats of other programs (SAM,
    // AIRR rearrangement, FASTA) before its warning. The deferred errors of `resolve` are
    // those of the composition-based statistics handler and of `Validate`, which come later.
    app::formatting_handler_check(&choice, Program::Blastx)?;
    // The warnings of the hit list size and of an ignored -num_descriptions
    // (`blastinput::report_counts`).
    stderr.write_all(&warnings_of(&args, choice.number))?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-175
    // ```c++
    //     catch (const blast::CInputException& e) {                               \
    //         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
    //         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
    //         exit_code = BLAST_INPUT_ERROR;                                      \
    // ```
    let mut options = resolution.map_err(|e| NativeError {
        exit: 1,
        message: format!(
            "BLAST query/options error: {e}\nPlease refer to the BLAST+ user manual.\n"
        ),
    })?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:195-200
    // ```c++
    //         /*** Initialize the database/subject ***/
    //         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
    //         CRef<CLocalDbAdapter> db_adapter;
    //         CRef<CScope> scope;
    //         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
    //                          db_adapter, scope);
    // ```
    // The database of `-db` opens after the arguments (the query and output files) and the
    // options are processed and before the query is read, as NCBI's `InitializeSubject`.
    let database = match (args.database_name(), &db_settings) {
        (Some(name), Some(settings)) => {
            let database =
                BlastxDatabase::open_restricted(name, settings, list.as_ref(), &mut stderr)?
                    .showing_gis(args.show_gis);
            crate::algorithm::common::database::note_archive_database_filter(
                archive.as_ref(),
                &args.db_filter,
                database.list_taxids(),
                database.filtered_empty(),
                args.searchsp.is_some_and(|space| space > 0),
                args.dbsize.is_some_and(|size| size > 0),
                "BLASTX",
            )?;
            Some(database)
        }
        _ => None,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1093-1104
    // ```c++
    // void CBlastTabularInfo::x_CheckTaxDB()
    // {
    //    	if( x_IsFieldRequested(eSubjectSciNames) ||
    //     ...
    //    		string resolved = SeqDB_ResolveDbPath("taxdb.bti");
    //    		if(resolved.empty())
    //    			ERR_POST(Warning << "Taxonomy name lookup from taxid requires installation of taxdb database with ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz");
    // ```
    options.filtered_empty = database
        .as_ref()
        .is_some_and(|database| database.filtered_empty());
    // A custom tabular field list that shows a plural name of the taxonomy database asks
    // whether `taxdb.bti` is on the search path: its formatters warn without it. A `-db` search
    // asks its database; a `-subject` search the search path of NCBI's databases (the
    // environment and the registry, read here for it alone), as BLASTP's and BLASTN's.
    if options.tabular.is_some_and(|spec| spec.checks_taxdb()) {
        options.taxdb_found = match &database {
            Some(database) => database.taxdb_found(),
            None => {
                let settings =
                    crate::blastinput::ncbi_environment::check_ncbi_application_settings("blastx")
                        .map_err(anyhow::Error::msg)?;
                crate::algorithm::common::database::subject_taxdb_found(&settings, "BLASTX")?
            }
        };
    }
    // LOSAT's refusals of a database for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let (Some(database), Some(serial)) = (
        database.as_ref(),
        app::ReportFormat::from_number(choice.number)
            .and_then(crate::report::serial::SerialFormat::of),
    ) {
        crate::report::serial::seqalign::check_database(
            serial,
            database.name(),
            &database.pairwise_database().definition,
            Program::Blastx,
        )?;
    }
    // LOSAT's refusal of a database name that BLAST XML (outfmt 5, 15 and 16) cannot write,
    // before the search (`blast_xml::check_database_name`).
    if let (Some(database), true) = (database.as_ref(), structured_xml) {
        super::blast_xml::check_database_name(database.name())?;
    }
    // LOSAT's other refusals of BLAST XML (a title with a byte NCBI's XML stream cannot write)
    // come before the first byte, though the query file is read batch by batch: the query
    // titles are checked by a read of the file before the search, and a run in which a
    // `-subject` title, the message of a query or a subject of a database with parsed Seq-ids
    // could still be refused later holds the document in memory until its epilog (`blast_xml::HeldOutput`; outfmt 13 and 14: the side files
    // too). Only with outfmt 5 and 13 to 16; the other formats read the file once.
    if let Some(held) = &held {
        // A query that is not a regular file (a pipe, `<(...)`, `/dev/stdin` from a pipe, a
        // FIFO) or `-query -` cannot be read a second time: the read before the search would
        // take every byte (DECISIONS 22). Refused before the search.
        if args.query.as_os_str() == "-"
            || !query_file.metadata().is_ok_and(|m| m.is_file())
            || !crate::blastinput::input_files::has_stream_position(&query_file)
        {
            bail!(
                "unsupported BLASTX query input for -outfmt {}: use a regular file",
                choice.number
            );
        }
        let query_file = File::open(&args.query).map_err(|_| inaccessible("query", &args.query))?;
        let messages_may_be_refused =
            super::blast_xml::prescan_queries(query_file, args.lcase_masking)?;
        // The documents of a `-db` search with a taxonomy database on the search path: a later
        // batch can still refuse the scientific name of a taxid (XML2 and JSON), or the damaged
        // database (`read_reported_taxids`, which needs the batch's hits; outfmt 5 looks the
        // batch's subjects up in its prefetch, `ReportedSubjects::prefetch_messages`, DBID item
        // T). The taxids of every sequence are looked up here, once, and the output is held only
        // when one of them would be refused (DECISIONS 32; outfmt 5 asks no scientific name); the
        // usual run, in which none is, writes batch by batch. Outfmt 5's prefetch is also refused
        // at the run's first prefetch where NCBI's taxonomy cache would repeat the offset error
        // (`BlastxDatabase::prefetch_refused`), after the prolog: the run holds then too.
        let prefetch = prefetches_subjects(Some(choice.number));
        let taxdb_may_refuse = database.as_ref().is_some_and(|database| {
            database.taxdb_found()
                && (database.taxdb_may_refuse(!prefetch) || prefetch && database.prefetch_refused())
        });
        // A database with parsed Seq-ids: a later batch can still refuse one of its subjects
        // (`ReportedSubjects::for_batch`), which only the batch's hits show.
        let parsed_may_be_refused = database
            .as_ref()
            .is_some_and(|database| database.parsed_seqids());
        if messages_may_be_refused
            || taxdb_may_refuse
            || parsed_may_be_refused
            || super::blast_xml::subjects_may_be_refused(&subjects)
        {
            held.hold();
        }
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-873
    // ```c++
    // 	char c;
    // 	CNcbiStreampos orig_p = in.tellg();
    // 	// Piped input
    // 	if(orig_p < 0)
    // 		return false;
    //
    // 	IOS_BASE::iostate orig_state = in.rdstate();
    // 	IOS_BASE::fmtflags orig_flags = in.setf(ios::skipws);
    //
    // 	if(! (in >> c))
    // 		return true;
    //
    // 	in.seekg(orig_p);
    // 	in.flags(orig_flags);
    // 	in.clear();
    // 	in.setstate(orig_state);
    //
    // 	return false;
    // ```
    // The prescan of the titles of outfmt 11 (below) reads the query file a second time: whether it
    // is a regular file is decided here, on the handle already open, before any second open.
    let query_is_regular = args.query.as_os_str() != "-" && archive::query_is_regular(&query_file);
    let mut query_stream = FastaStream::from_file(query_file);
    query_stream.parse_deflines = args.parse_deflines;
    query_stream.strand = crate::report::archive::bioseq::SeqLocStrand::of_option(&args.strand);
    if stream_is_empty(&mut query_stream) {
        stderr.write_all(b"Warning: [blastx] Query is Empty!\n")?;
        // NCBI returns before the formatter exists: no archive carries the message.
        if let Some(handle) = &archive {
            handle.take_messages();
        }
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:217-218
    // ```c++
    //         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
    //         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:85-91
    // ```c++
    //     // used for experimentation purposes
    //     char* batch_sz_str = getenv("BATCH_SIZE");
    //     if (batch_sz_str) {
    //         retval = NStr::StringToInt(batch_sz_str);
    //         _TRACE("DEBUG: Using query batch size " << retval);
    //         return retval;
    //     }
    // ```
    // The letters of a reader / report batch: `BATCH_SIZE` when it is set, else 10002. The
    // search itself keeps units of at most 10002 letters (`BATCH_SIZE` the constant), however
    // large a reader batch is: a reader batch is split into units that are searched as the
    // batches of a run without `BATCH_SIZE` are, and its report is written when its last unit
    // is searched.
    let batch_size =
        crate::blastinput::app::query_batch_size("BLASTX", BATCH_SIZE as u32)? as usize;
    // LOSAT's refusal of a query title the archive cannot write comes before its first byte,
    // though the query file is read batch by batch: the titles are checked by a read of the
    // file before the search.
    if archive.is_some() {
        if !query_is_regular {
            archive::refuse_irregular_query()?;
        }
        archive::prescan_query_titles(&args.query, args.lcase_masking, batch_size)?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
    // ```c++
    //         catch(CBlastException & e ) {
    //         	// Skip bad subject sequence
    //         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
    //         		seqblk_vec->push_back(subj);
    //         		string warning = kEmptyStr;
    //         		const CSeq_id *  id = subjects.GetSeqId(i);
    //         		string title = subjects.GetTitle(i);
    //         		if(id != NULL) {
    //         			warning = id->GetSeqIdString() + " ";
    //         		}
    //         		warning += subjects.GetTitle(i);
    //         		if(warning != kEmptyStr){
    //         			warning += ": ";
    //         		}
    //         		warning += "Subject sequence contains no data";
    //         		ERR_POST(Warning << warning);
    //         		continue;
    // ```
    for s in subjects.iter().filter(|s| s.sequence.is_empty()) {
        writeln!(
            stderr,
            "Warning: [blastx] {} {}: Subject sequence contains no data",
            String::from_utf8_lossy(s.warning_id()),
            s.title
        )?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:251-254
    // ```c++
    // 	if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
    //         	formatter.SetBaseFile(args[kArgOutput].AsString());
    //         }
    //         formatter.PrintProlog();
    // ```
    // The XInclude formats written to standard output fail after the formatter is made (a
    // format 13 or 14 to a file is LOSAT's rejection, above).
    app::xinclude_check_for(&choice, args.out.as_deref())?;
    // LOSAT_X_DBIDXREAD acts on a database's index files only; the new side is marked where a
    // volume is opened (`seqdb::Volume::open`).
    if database.is_none() && crate::seqdb::x_dbidxread() != crate::seqdb::XDbIdxRead::Off {
        crate::utils::xreach::x_reach!("LOSAT_X_DBIDXREAD", false);
    }
    let pairwise_database = database.as_ref().map(BlastxDatabase::pairwise_database);
    let no_reported = ReportedSubjects::default();
    let records_parsed = if args.parse_deflines {
        report::RecordsParsed::of(&subjects, args.show_gis, options.outfmt == 0)?
    } else {
        report::RecordsParsed::default()
    };
    let shown_run = match (&database, &pairwise_database) {
        (Some(database), Some(pairwise)) => ShownSubjects::Database {
            name: database.name(),
            pairwise,
            reported: &no_reported,
        },
        _ => ShownSubjects::Fasta {
            path: subject_path.expect("a subject file without -db"),
            parsed: &records_parsed,
        },
    };
    super::x_timing::run(XRun::Startup, x_t_run);
    let x_t = super::x_timing::lap();
    // The writer of the XML2 or JSON document (outfmt 15 and 16), which the run keeps across the
    // query batches; its start is written with the prolog, as NCBI writes `PrintProlog`'s.
    let xml2 = std::cell::RefCell::new(super::blast_xml::xml2_kind(options.outfmt).map(|kind| {
        let files: RunSideFiles = match &side_files {
            Some(side_files) => Box::new(side_files.clone()),
            None => Box::new(NoSideFiles),
        };
        super::blast_xml::xml2_writer(kind, &options, &subjects, &shown_run, files)
    }));
    report::write_prolog(&mut *writer, &mut *stderr, &options, &subjects, &shown_run)
        .map_err(report_error)?;
    if let Some(xml2) = xml2.borrow_mut().as_mut() {
        let mut start = Vec::new();
        xml2.start(&mut start);
        writer.write_all(&start).map_err(output_error)?;
    }
    // NCBI writes the prolog out before it searches the first batch, so the warnings of the
    // batch come after it when both streams go to one place.
    writer.flush().map_err(output_error)?;
    super::x_timing::run(XRun::Prolog, x_t);
    let x_t = super::x_timing::lap();
    // The database's sequences, checked once here and read in blocks where they are searched
    // (`x_db_stream.rs`); an error belongs to the first batch with a valid context
    // (`Subjects::Unreadable`, `search_error`). No NCBI counterpart for the blocks: NCBI maps the
    // volume files.
    let (x_longest, unreadable) = match &database {
        None => (None, None),
        Some(database) => match database.x_check_stores() {
            Ok(longest) => (Some(longest), None),
            Err(error) => (None, Some(error)),
        },
    };
    let x_stream = match (&database, x_longest) {
        (Some(database), Some(longest)) => Some(super::x_db_stream::XStreamDb::new(
            database.x_seqdb(),
            longest,
        )),
        _ => None,
    };
    super::x_timing::run(XRun::ReadStores, x_t);
    let unreadable = std::cell::RefCell::new(unreadable);
    let search_subjects = match (&database, &x_stream) {
        (None, _) => Subjects::Fasta(&subjects),
        (Some(_), Some(stream)) => Subjects::Stream(super::x_db_stream::XStream {
            db: stream,
            block: None,
        }),
        (Some(_), None) => Subjects::Unreadable,
    };
    let search_error = |error: anyhow::Error| -> anyhow::Error {
        // LOSAT's refusal, not an error of NCBI's engine.
        if error.is::<super::parameters::SearchSpaceOverflow>() {
            return error;
        }
        if error.is::<UnreadableSubjects>() {
            let kept = unreadable
                .borrow_mut()
                .take()
                .expect("the database's error, raised once");
            return database_error(kept);
        }
        // A read of the database that failed during the search,
        // reported as the check's error for the same read.
        match error.downcast::<super::x_db_stream::XReadError>() {
            Ok(super::x_db_stream::XReadError(kept)) => database_error(kept),
            Err(error) => engine(error),
        }
    };
    // The reports of one query batch: with `-db`, the titles of its reported sequences are
    // read first and its hit lists renamed to ranks (`ReportedSubjects::for_batch`); NCBI's
    // reader error there stops the run at this batch, with no epilog.
    // A refusal of a title (BLAST XML) found while a batch is reported drops what the run holds.
    let refused = |error: anyhow::Error| -> anyhow::Error {
        let error = report_error(error);
        if error.is::<super::blast_xml::XmlRefusal>() {
            if let Some(held) = &held {
                held.discard();
            }
        }
        error
    };
    let render_batch = |writer: &mut BufWriter<Box<dyn Write>>,
                        diagnostics: &mut dyn Write,
                        queries: &[FastaRecord],
                        results: &mut Vec<BatchResults>,
                        queries_before: usize|
     -> Result<()> {
        if let Some(handle) = &archive {
            // The BLAST archive (outfmt 11): one document per query batch, written from the
            // data of outfmt 8's Seq-aligns.
            for batch in results.iter_mut() {
                for lists in &mut batch.queries {
                    for list in lists {
                        list.sort_evalue();
                    }
                }
            }
            let reported;
            let (subject_records, subject_names, call_subjects): (_, _, &[FastaRecord]) =
                match (&database, &pairwise_database) {
                    (Some(database), Some(pairwise)) => {
                        reported = match search_subjects {
                            Subjects::Unreadable => ReportedSubjects::default(),
                            search_subjects => ReportedSubjects::for_batch(
                                database,
                                search_subjects,
                                results,
                                options.gapped,
                                false,
                                false,
                                false,
                                false,
                            )
                            .map_err(|error| match error.downcast::<HeaderReadError>() {
                                Ok(header) => header.error,
                                Err(error) => {
                                    match error.downcast::<super::x_db_stream::XReadError>() {
                                        Ok(super::x_db_stream::XReadError(kept)) => {
                                            database_error(kept)
                                        }
                                        Err(error) => error,
                                    }
                                }
                            })?,
                        };
                        (
                            None,
                            super::seqalign::SeqalignSubjects::Database {
                                reported: &reported,
                                title: &pairwise.definition,
                            },
                            reported.records(),
                        )
                    }
                    _ => (
                        Some(&subjects[..]),
                        super::seqalign::SeqalignSubjects::Records,
                        &subjects[..],
                    ),
                };
            return archive::write_blastx_archive(
                &archive::ArchiveCall {
                    handle,
                    options: &options,
                    queries,
                    subjects: call_subjects,
                    names: subject_names,
                    subject_records,
                    batches: results,
                },
                writer,
            )
            .map_err(archive_error);
        }
        match (&database, &pairwise_database) {
            (Some(database), Some(pairwise)) => {
                let reported = match search_subjects {
                    Subjects::Unreadable => ReportedSubjects::default(),
                    search_subjects => {
                        let x_t = super::x_timing::lap();
                        let reported = ReportedSubjects::for_batch(
                            database,
                            search_subjects,
                            results,
                            options.gapped,
                            matches!(options.outfmt, 6 | 7 | 10 | 20),
                            options.outfmt == 0,
                            matches!(options.outfmt, 13..=16),
                            matches!(options.outfmt, 0 | 5 | 13..=16),
                        );
                        super::x_timing::run(XRun::Titles, x_t);
                        let reported = reported
                            .map_err(|error| match error.downcast::<HeaderReadError>() {
                                Ok(header) => header.error,
                                // A failed read of a record's residues.
                                Err(error) => {
                                    match error.downcast::<super::x_db_stream::XReadError>() {
                                        Ok(super::x_db_stream::XReadError(kept)) => {
                                            database_error(kept)
                                        }
                                        Err(error) => error,
                                    }
                                }
                            })
                            .inspect_err(|error| {
                                // A refusal of the XML2 document, or of a subject with parsed
                                // Seq-ids, drops what the run holds.
                                if error.is::<super::blast_xml::XmlRefusal>()
                                    || error
                                        .is::<crate::algorithm::common::database::ParsedRefusal>()
                                {
                                    if let Some(held) = &held {
                                        held.discard();
                                    }
                                }
                            })?;
                        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:286-291
                        // ```c++
                        //                 BlastFormatter_PreFetchSequenceData(*results, scope,
                        //                 		                            fmt_args->GetFormattedOutputChoice());
                        //             	ITERATE(CSearchResultSet, result, *results) {
                        //                	    formatter.PrintOneResultSet(**result, query_batch);
                        //             	}
                        // ```
                        // Before the reports of a batch NCBI fetches every subject that the batch
                        // reports, and with `taxdb.bti` on the search path the taxonomy database
                        // writes its messages there (`ReportedSubjects::prefetch_messages`: the
                        // texts of opening it at the run's first lookup, a volume's lookup of a
                        // taxid once per run): after the reader's messages of the reader batch,
                        // which every batch loop writes before it calls this, and before the
                        // warnings of its first query (`report::render_queries`). Outfmt 0, 5, 6
                        // and 7 (`database::prefetches_subjects`); XML2 and JSON write theirs at
                        // each `Hit` (`TaxdbPrinter`). Without `taxdb.bti` nothing is read. A
                        // refusal comes after the reports of the batches before, as the refusals
                        // of the batch's titles above (line DBID, DECISIONS 9 (9)); a refusal of
                        // the taxonomy database (`TaxonomyRefusal`) drops what an outfmt 5 run
                        // holds, and the run holds where one could come (`taxdb_may_refuse`
                        // above), as XML2's does.
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
                        // Outfmt 10 and 20 read the subjects as 6 and 7 do; the rows' own
                        // lookups of the taxonomy names come with the rows (`tabular`).
                        let prefetch_format = match options.outfmt {
                            10 | 20 => 6,
                            outfmt => i32::from(outfmt),
                        };
                        if prefetches_subjects(Some(prefetch_format)) && database.taxdb_found() {
                            let messages = reported
                                .prefetch_messages(database, results, options.gapped)
                                .inspect_err(|error| {
                                    if error.is::<TaxonomyRefusal>() {
                                        if let Some(held) = &held {
                                            held.discard();
                                        }
                                    }
                                })?;
                            diagnostics.write_all(&messages)?;
                        }
                        // The values of a custom tabular field list that reads more of the
                        // definition lines than the default report (`sall*`) or the taxonomy
                        // (`staxid`, `staxids`, the names), read before the batch's first byte so
                        // that their refusals come before it, after the prefetch (which opens the
                        // taxonomy database first, as NCBI's does).
                        let mut reported = reported;
                        if let Some(spec) = options.tabular {
                            reported.read_tabular(database, spec, results, options.gapped)?;
                        }
                        reported
                    }
                };
                report::render_queries(
                    writer,
                    diagnostics,
                    &options,
                    queries,
                    reported.records(),
                    &ShownSubjects::Database {
                        name: database.name(),
                        pairwise,
                        reported: &reported,
                    },
                    results,
                    queries_before,
                    xml2.borrow_mut().as_mut(),
                )
                .map_err(refused)
            }
            _ => report::render_queries(
                writer,
                diagnostics,
                &options,
                queries,
                &subjects,
                &shown_run,
                results,
                queries_before,
                xml2.borrow_mut().as_mut(),
            )
            .map_err(refused),
        }
    };
    // NCBI: a batch size of 0 makes the first batch empty, and its query factory refuses it.
    if batch_size == 0 && !query_stream.at_eof() {
        return Err(engine("Empty CBlastQueryVector"));
    }
    let mut processed = 0usize;
    let x_t_loop = super::x_timing::lap();
    {
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:260-261
        // ```c++
        // 	    BLAST_PROF_START( APP.LOOP.PRE );
        //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
        // ```
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:289-290
        // ```c++
        //             	ITERATE(CSearchResultSet, result, *results) {
        //                	    formatter.PrintOneResultSet(**result, query_batch);
        // ```
        // The reader and formatter use the same diagnostic stream at distinct
        // times. RefCell only serializes Rust borrows; it does not defer events.
        let diagnostics = std::cell::RefCell::new(&mut stderr);
        // The batch loop as a function of the pool its searches run on, so that one pool can serve
        // all batches.
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-291
        // ```c
        //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
        // ...
        //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
        // ...
        //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
        //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
        //                 results = lcl_blast.Run();
        // ...
        //             	ITERATE(CSearchResultSet, result, *results) {
        //                	    formatter.PrintOneResultSet(**result, query_batch);
        //             	}
        // ```
        // NCBI runs this loop: take the next batch of queries, search it, print its results.
        // x_run_batches is the same loop (the code below the batched branch is the original loop
        // body) as a closure over the pool its searches run on, so that one pool can serve all
        // batches.
        let mut x_run_batches = |x_pool: Option<&crate::utils::threading::SearchPool<'_>>| {
            // With a pool of 2+ threads: search several query batches at a time.
            //
            // A batch (`BATCH_SIZE` bases of queries) is searched on its own:
            // nothing is carried from one batch to the next except the order of
            // what is written. A file of many short queries is therefore a long
            // row of small searches whose serial parts (query masking, lookup
            // table, hand-over between stages) leave the other threads idle. Here
            // the reader's events are queued instead of acted on: a warning as the
            // bytes it would have written, a batch as a copy of its records. When
            // enough batches are queued they are searched side by side on the
            // pool, and the queue is then replayed in order: warnings are written,
            // each batch is reported exactly where the unbatched loop reported it.
            // A batch that fails ends the replay there, as it ended the loop; what
            // the reader delivered before it failed itself is replayed first.
            // A batch with a long query is split into chunks that use the whole
            // pool already; it is flushed at once, never queued with others.
            // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:345-357
            // ```c
            // 		while (master_node.Processing()) {
            // 			if (!input.AtEOF()) {
            // 			 	if (!master_node.IsFull()) {
            // ...
            // 					int num_q = input.GetQueryBatch(qb, q_index);
            // 					if (num_q > 0) {
            // 						CBlastNodeMailbox * mb(new CBlastNodeMailbox(chunk_num, master_node.GetBuzzer()));
            // 						CBlastxNode * t(new CBlastxNode(chunk_num, GetArguments(), args, m_Bah, qb, q_index, num_q, mb));
            // 						master_node.RegisterNode(t, mb);
            // 						chunk_num ++;
            // ```
            // NCBI has a similar mode (-mt_mode 1, x_RunMTBySplitQuery): several query batches are
            // searched at the same time as nodes of a master node. Here the batches are queued and
            // replayed in input order; the results printed are those of the one-batch-at-a-time
            // loop above.
            #[cfg(feature = "parallel")]
            {
                // Without a pool of 2+ threads the batches are searched one at a time;
                // LOSAT_X_BXFLOW acts only with the queued batches.
                let x_pooled = x_pool.is_some_and(|pool| pool.enabled());
                if !x_pooled {
                    crate::utils::xreach::x_reach!("LOSAT_X_BXBATCH", false);
                }
                if super::runtime::x_bx_flow() {
                    crate::utils::xreach::x_reach!("LOSAT_X_BXFLOW", x_pooled);
                }
            }
            #[cfg(feature = "parallel")]
            // With `-parse_deflines` a batch is searched and reported before the next is read: a
            // reported database sequence is refused by the ids of the queries read so far
            // (`common::database::read_header`), which must not depend on the threads.
            if let Some(pool) = x_pool.filter(|pool| pool.enabled() && !args.parse_deflines) {
                use rayon::prelude::*;
                enum XEvent {
                    Warn(Vec<u8>),
                    Batch(XUnit),
                }
                let x_events = std::cell::RefCell::new(Vec::<XEvent>::new());
                let x_failed = std::cell::Cell::new(false);
                let mut pending = ReaderBatchResults::default();
                let x_wave = pool.threads() * 4;
                // A unit of more than two search units' letters ends in a query that is longer than
                // one: it is searched with the pool inside, never queued with others. The query
                // is split into chunks that use the pool when it is the only query of its reader
                // batch; in a reader batch of several queries (NCBI splits none) its subjects use
                // the pool. Measured (2 queries of 280 kb, `BATCH_SIZE=1000000`, 4 and 35
                // subjects, 8 and 16 threads): 0.10-0.16 s, against 0.17-0.26 s for the same
                // queries split.
                let x_large = |queries: &[super::input::FastaRecord]| {
                    queries.iter().map(|q| q.sequence.len()).sum::<usize>() > 2 * BATCH_SIZE
                };
                // A reader batch becomes its units; a batch of several is prepared first (phase
                // 1, `reader_batch.rs`): by the waves at once, on the pool, by the flow as tasks.
                let x_split = |queries: &[super::input::FastaRecord]| {
                    let units = XUnit::split(queries, &options, false);
                    XUnit::prepare_now(&units, &|batch, ranges| {
                        prepare_reader_batch(batch, ranges, &options, search_subjects, pool)
                    });
                    units
                };
                // EXPERIMENT (LOSAT_X_BXFLOW): the queued batches as a flow instead of waves.
                //
                // The waves below search `x_wave` batches and then wait for the last of them
                // before any is written or the next wave is read; at the end of the input the
                // remaining batches are searched one after another. Here every batch becomes a
                // task of its own once the reader has delivered as many batches after it as there
                // are threads, and this thread (worker 0 of the pool) writes the queued warnings
                // and reports in the reader's order as the results come in. The last batches of
                // the input (at most one per thread, the waves' tail) are each searched with the
                // pool inside, as the waves search their tail, but side by side. At most `x_wave`
                // batches are read and not yet written; at that bound, and at the end of the
                // input, this thread runs pending tasks itself and otherwise waits for a result.
                // Each batch is searched by a function the waves use, so the reports, the
                // warnings, their order, and the stop at a batch that fails are those of the
                // waves; only the time at which a batch is searched changes.
                // No NCBI counterpart: NCBI's -mt_mode 1 (blastx_app.cpp x_RunMTBySplitQuery, cited
                // above x_flush) also searches several query batches at once and prints them in
                // order; the flow changes only LOSAT's scheduling of the searches above.
                if super::runtime::x_bx_flow() {
                    // Dispatch point (LOSAT_X_BXUNITCAP=K, in the strict set): each batch searched
                    // with the pool inside (the subject flow) instead of one batch per thread, and
                    // at most K at a time: the flow spawns a batch once one batch follows it and
                    // while fewer than K are searched, and reads on while fewer than
                    // `flow_bound(K)` are unwritten (blastx/x_bx_mem.rs).
                    let x_cap =
                        super::search::x_bx_mem::flow_unit_cap(pool, &options, search_subjects);
                    let (x_bound, x_lookahead, x_max_live) = match x_cap {
                        Some(k) => (super::search::x_bx_mem::flow_bound(k), 1, k),
                        None => (x_wave, pool.threads(), usize::MAX),
                    };
                    let x_ahead =
                        x_cap.is_some() && super::search::x_bx_mem::prepares_batches_ahead();
                    return x_flow(
                        x_bound,
                        x_lookahead,
                        x_max_live,
                        &mut query_stream,
                        args.lcase_masking,
                        batch_size,
                        archive.is_some(),
                        // Dispatch point (LOSAT_X_BXCAPFILL, in the strict set): in the capped
                        // flow a reader batch of one unit is prepared by a task of its own before
                        // its search, outside the K searched at a time (blastx/x_bx_mem.rs).
                        &|queries: &[super::input::FastaRecord]| {
                            XUnit::split(queries, &options, x_ahead)
                        },
                        &|queries: &[super::input::FastaRecord],
                          layout: Option<&std::sync::Arc<ReaderBatchShape>>| match layout {
                            Some(layout) => prepare_unit(queries, layout, &options, search_subjects),
                            None => prepare_lone_unit(queries, &options),
                        },
                        &|unit: &XUnit, alone: bool| {
                            let queries = unit.queries();
                            let large = x_cap.is_some() || alone || x_large(queries);
                            crate::utils::xreach::x_reach!("LOSAT_X_BXBATCH", !large);
                            if large {
                                let _x_live = x_cap
                                    .map(|_| super::search::x_bx_mem::XCappedBatch::enter(pool));
                                super::runtime::x_search_unit_in_pool(
                                    queries,
                                    search_subjects,
                                    &options,
                                    pool,
                                    unit.take_prep(),
                                )
                            } else {
                                super::runtime::x_search_unit_serial(
                                    queries,
                                    search_subjects,
                                    &options,
                                    unit.take_prep(),
                                )
                            }
                        },
                        &mut |reader_batch: &[super::input::FastaRecord],
                              start: usize,
                              last: bool,
                              results: Result<Vec<BatchResults>>|
                         -> Result<()> {
                            let results = results.map_err(search_error)?;
                            if let Some(mut results) = pending.add(start, last, results) {
                                render_batch(
                                    &mut *writer,
                                    &mut **diagnostics.borrow_mut(),
                                    reader_batch,
                                    &mut results,
                                    processed,
                                )?;
                                writer.flush().map_err(output_error)?;
                                processed += reader_batch.len();
                            }
                            Ok(())
                        },
                        &mut |bytes: &[u8]| -> Result<()> {
                            diagnostics.borrow_mut().write_all(bytes)?;
                            Ok(())
                        },
                    );
                }
                // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-291
                // ```c
                //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
                // ...
                //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
                // ...
                //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
                //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
                //                 results = lcl_blast.Run();
                // ...
                //             	ITERATE(CSearchResultSet, result, *results) {
                //                	    formatter.PrintOneResultSet(**result, query_batch);
                //             	}
                // ```
                // x_flush searches the queued batches (side by side, or one after another with the
                // whole pool inside each) and then writes warnings and results in the order the
                // reader produced them, which is the order of the loop above.
                let mut x_flush = |events: Vec<XEvent>| -> Result<()> {
                    let batches: Vec<&XUnit> = events
                        .iter()
                        .filter_map(|event| match event {
                            XEvent::Batch(unit) => Some(unit),
                            XEvent::Warn(_) => None,
                        })
                        .collect();
                    // Fewer batches than threads (the end of the input): each one
                    // uses the pool inside, as the unbatched loop does.
                    let x_t = super::x_timing::lap();
                    let x_tail = batches.len() < pool.threads();
                    crate::utils::xreach::x_reach!("LOSAT_X_BXBATCH", !x_tail);
                    let searched = if x_tail {
                        batches
                            .iter()
                            .map(|unit| {
                                super::runtime::x_search_unit_in_pool(
                                    unit.queries(),
                                    search_subjects,
                                    &options,
                                    pool,
                                    unit.take_prep(),
                                )
                            })
                            .collect::<Vec<_>>()
                    } else {
                        pool.install(|| {
                            batches
                                .par_iter()
                                .with_max_len(1)
                                .map(|unit| {
                                    let x_t = super::x_timing::lap();
                                    // A long unit uses the pool inside (a single long query
                                    // is split into chunks, otherwise its subjects use the
                                    // pool); a short one is one serial job.
                                    let queries = unit.queries();
                                    let searched = if x_large(queries) {
                                        super::runtime::x_search_unit_in_pool(
                                            queries,
                                            search_subjects,
                                            &options,
                                            pool,
                                            unit.take_prep(),
                                        )
                                    } else {
                                        super::runtime::x_search_unit_serial(
                                            queries,
                                            search_subjects,
                                            &options,
                                            unit.take_prep(),
                                        )
                                    };
                                    super::x_timing::run(XRun::WaveBusy, x_t);
                                    searched
                                })
                                .collect::<Vec<_>>()
                        })
                    };
                    // Diagnostic (LOSAT_TIMING): the wall time of the searches. No NCBI
                    // counterpart.
                    let x_wall = super::x_timing::run(XRun::Search, x_t);
                    let (x_count, x_wall_stage) = if x_tail {
                        (XRun::TailBatches, XRun::TailWall)
                    } else {
                        super::x_timing::run_add(XRun::Waves, 1);
                        (XRun::WaveBatches, XRun::WaveWall)
                    };
                    super::x_timing::run_add(x_count, batches.len() as u64);
                    super::x_timing::run_add(x_wall_stage, x_wall);
                    let x_t = super::x_timing::lap();
                    let mut searched = searched.into_iter();
                    for event in &events {
                        match event {
                            XEvent::Warn(bytes) => diagnostics.borrow_mut().write_all(bytes)?,
                            XEvent::Batch(unit) => {
                                let results = searched
                                    .next()
                                    .expect("one search per queued batch")
                                    .map_err(search_error)?;
                                if let Some(mut results) =
                                    pending.add(unit.range.start, unit.last, results)
                                {
                                    render_batch(
                                        &mut *writer,
                                        &mut **diagnostics.borrow_mut(),
                                        &unit.batch,
                                        &mut results,
                                        processed,
                                    )?;
                                    writer.flush().map_err(output_error)?;
                                    processed += unit.batch.len();
                                }
                            }
                        }
                    }
                    super::x_timing::run(XRun::Render, x_t);
                    Ok(())
                };
                let parsed = parse_fasta_stream_keeping(
                    &mut query_stream,
                    false,
                    args.lcase_masking,
                    batch_size,
                    archive.is_some(),
                    &mut |queries| {
                        if queries.is_empty() {
                            return Err(engine("Empty CBlastQueryVector"));
                        }
                        // The reader batch is queued as its units (at most 10002 letters
                        // each), which are searched as the batches of a run without
                        // `BATCH_SIZE` are; its report is written at its last unit.
                        for unit in x_split(queries) {
                            let large = x_large(unit.queries());
                            let events = {
                                let mut events = x_events.borrow_mut();
                                events.push(XEvent::Batch(unit));
                                let queued = events
                                    .iter()
                                    .filter(|event| matches!(event, XEvent::Batch(_)))
                                    .count();
                                if !large && queued < x_wave {
                                    continue;
                                }
                                std::mem::take(&mut *events)
                            };
                            let flushed = x_flush(events);
                            if flushed.is_err() {
                                x_failed.set(true);
                                return flushed;
                            }
                        }
                        Ok(())
                    },
                    &mut |warning| {
                        let mut bytes = Vec::new();
                        reader_warning(warning, &mut bytes)?;
                        x_events.borrow_mut().push(XEvent::Warn(bytes));
                        Ok(())
                    },
                );
                if !x_failed.get() {
                    x_flush(x_events.take())?;
                }
                return parsed;
            }
            parse_fasta_stream_keeping(
                &mut query_stream,
                false,
                args.lcase_masking,
                batch_size,
                archive.is_some(),
                &mut |queries| {
                    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-262
                    // ```c++
                    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
                    // 	    BLAST_PROF_START( APP.LOOP.PRE );
                    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
                    //             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
                    // ```
                    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
                    // ```c++
                    // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
                    //     : m_QueryVector(& queries)
                    // {
                    //     if (queries.Empty()) {
                    //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
                    //     }
                    // ```
                    if queries.is_empty() {
                        return Err(engine("Empty CBlastQueryVector"));
                    }
                    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
                    // ```c
                    //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
                    //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
                    //                 results = lcl_blast.Run();
                    // ```
                    // Dispatch point: with a kept pool the batch is searched on it; without one,
                    // search_internal runs as before. Both searches compute the same result (the
                    // pool only changes who runs the work).
                    let x_t = super::x_timing::lap();
                    let searched = match x_pool {
                        Some(pool) => super::runtime::x_search_internal_in_pool(
                            queries,
                            search_subjects,
                            &options,
                            pool,
                        ),
                        None => search_internal(queries, search_subjects, &options),
                    };
                    super::x_timing::run(XRun::Search, x_t);
                    let x_t = super::x_timing::lap();
                    let mut results = searched.map_err(search_error)?;
                    reader_batch_skipped(&mut results);
                    render_batch(
                        &mut *writer,
                        &mut **diagnostics.borrow_mut(),
                        queries,
                        &mut results,
                        processed,
                    )?;
                    writer.flush().map_err(output_error)?;
                    super::x_timing::run(XRun::Render, x_t);
                    processed += queries.len();
                    Ok(())
                },
                &mut |warning| reader_warning(warning, &mut **diagnostics.borrow_mut()),
            )
        };
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
        // ```c
        //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
        //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
        //                 results = lcl_blast.Run();
        // ```
        // Dispatch point: NCBI builds a CLocalBlast, and so its worker threads, for every batch.
        // With 2+ threads one pool is built for the whole run and every batch is searched on
        // it. Reuse of threads only.
        crate::utils::xreach::x_reach!("LOSAT_X_BXPOOL", options.num_threads > 1);
        let x_batches = if options.num_threads > 1 {
            // A pool that cannot be built is the engine error it was when each
            // batch built its own.
            match crate::utils::threading::x_with_search_pool_local(
                options.num_threads as usize,
                "BLASTX",
                |pool| {
                    // The run's counts of LOSAT_X_BXUNITCAP (scheduling only), on its pool.
                    #[cfg(feature = "parallel")]
                    let x_caps = super::search::x_bx_mem::XCapRun::default();
                    #[cfg(feature = "parallel")]
                    let pool = &x_caps.attach(pool);
                    Ok(x_run_batches(Some(pool)))
                },
            ) {
                Ok(batches) => batches,
                Err(error) => Err(engine(error)),
            }
        } else {
            x_run_batches(None)
        };
        x_batches.map_err(input_error)?;
    }
    super::x_timing::run(XRun::Loop, x_t_loop);
    let x_t = super::x_timing::lap();
    report::write_epilog(
        &mut *writer,
        &mut *stderr,
        &options,
        &subjects,
        &shown_run,
        processed,
        xml2.borrow_mut().as_mut(),
    )
    .map_err(report_error)?;
    writer.flush().map_err(output_error)?;
    // The document a hold kept: written at the end of the run, with the epilog.
    if let Some(held) = &held {
        held.release().map_err(report_error)?;
    }
    super::x_timing::run(XRun::Epilog, x_t);
    super::x_timing::run(XRun::Total, x_t_run);
    // Diagnostic (LOSAT_TIMING): the run's share of the common line (x_timing.rs).
    super::x_timing::common_run();
    super::x_timing::write(&mut *stderr, options.num_threads as usize)?;
    Ok(())
}

/// Runs the search; NCBI writes an error archive after the archives of the batches before when
/// the application stops with one of its own errors, and when the run ends with messages that no
/// archive carried (`PrintErrorArchive` after its `CATCH_ALL`, whatever the status), into the
/// output the run opened.
pub fn run(args: &BlastxArgs) -> Result<()> {
    let mut output: Option<BufWriter<Box<dyn Write>>> = None;
    let mut archive: Option<crate::report::archive::ArchiveHandle> = None;
    let result = run_inner(args, &mut output, &mut archive);
    if let (Some(handle), Some(writer)) = (&archive, output.as_mut()) {
        match &result {
            Err(error) => {
                if let Some(lines) = archive::error_lines(error, false) {
                    let _ = crate::report::archive::write_error_archive(handle, &lines, writer);
                }
            }
            Ok(()) => {
                if handle.has_messages() {
                    let _ = crate::report::archive::write_error_archive(handle, &[], writer);
                }
            }
        }
        let _ = writer.flush();
    }
    result
}

/// The error archive of a run that stopped in `Run` before its output was opened (the taxid or
/// seqid list, the files of `-subject` and `-query`, the subjects): into `-out`, or standard output; nothing when that
/// cannot be opened (NCBI's call swallows that).
fn write_run_error_archive(
    args: &BlastxArgs,
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
    let Some(lines) = archive::error_lines(error, true) else {
        return;
    };
    let destination: Option<Box<dyn Write>> =
        match args.out.as_deref().filter(|p| *p != Path::new("-")) {
            Some(path) => File::create(path)
                .ok()
                .map(|file| Box::new(file) as Box<dyn Write>),
            None => Some(Box::new(io::stdout().lock())),
        };
    if let Some(mut destination) = destination {
        let _ = crate::report::archive::write_error_archive(handle, &lines, &mut destination);
        let _ = destination.flush();
    }
}
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:989-1013
// ```c++
//         FASTA_WARNING_EX(LineNumber(),
//             "CFastaReader: Hyphens are invalid and will be ignored around line " << LineNumber(),
//             ILineError::eProblem_IgnoredResidue,
//             kEmptyStr, kEmptyStr, "-" );
//     }
//
//     // before throwing, be sure that we're in a valid state so that callers can
//     // parse multiple lines and get the invalid residues in all of them.
//
//     if( ! bad_pos_vec.empty() ) {
//         if (TestFlag(fValidate)) {
//                         NCBI_THROW2(CBadResiduesException, eBadResidues,
//                 "CFastaReader: There are invalid " + x_NucOrProt() + "residue(s) in input sequence",
//                 CBadResiduesException::SBadResiduePositions( m_BestID, bad_pos_vec, bad_pos_line_num ) );
//         } else {
//             stringstream warn_strm;
//             warn_strm << "FASTA-Reader: Ignoring invalid " << x_NucOrProt()
//                 << "residues at position(s): ";
//             CBadResiduesException::SBadResiduePositions(
//                 m_BestID, bad_pos_vec, bad_pos_line_num ).ConvertBadIndexesToString(warn_strm);
//
//             FASTA_WARNING(0,
//                 warn_strm.str(),
//                 ILineError::eProblem_InvalidResidue,
//                 kEmptyStr );
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2196-2217
// ```c++
// void CFastaReader::PostWarning(
//             ILineErrorListener * pMessageListener,
//             EDiagSev _eSeverity, size_t _uLineNum, CTempString _MessageStrmOps, CObjReaderParseException::EErrCode _eErrCode, ILineError::EProblem _eProblem, CTempString _sFeature, CTempString _sQualName, CTempString _sQualValue) const
// {
//     if (find(m_ignorable.begin(), m_ignorable.end(), _eProblem) != m_ignorable.end())
//         // this is a problem that should be ignored
//         return;
//
//     string sSeqId = ( m_BestID ? m_BestID->AsFastaString() : kEmptyStr);
//     AutoPtr<CObjReaderLineException> pLineExpt(
//         CObjReaderLineException::Create(
//         (_eSeverity), static_cast<unsigned>(_uLineNum),
//         _MessageStrmOps,
//         (_eProblem),
//         sSeqId, (_sFeature),
//         (_sQualName), (_sQualValue),
//         _eErrCode) );
//     if ( ! pMessageListener && (_eSeverity) <= eDiag_Warning ) {
//         LOG_POST_X(1, Warning << pLineExpt->Message());
//     } else if ( ! pMessageListener || ! pMessageListener->PutError( *pLineExpt ) )
//     {
//         throw CObjReaderParseException(DIAG_COMPILE_INFO, 0, _eErrCode, _MessageStrmOps, _uLineNum, _eSeverity);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_exception.cpp:106-137
// ```c++
//                 ++iRangesFound;
//                 continue;
//             }
//
//             const TSeqPos last_idx = rangesFound.back().second;
//             if( idx == (last_idx+1) ) {
//                 // extend previous range
//                 ++rangesFound.back().second;
//                 continue;
//             }
//
//             if( iRangesFound >= maxRanges ) {
//                 break;
//             }
//
//             // create new range
//             rangesFound.push_back(TRange(idx, idx));
//             ++iRangesFound;
//         }
//
//         // turn the ranges found on this line into a string
//         out << line_prefix << "On line " << lineNum << ": ";
//         line_prefix = ", ";
//
//         const char *pos_prefix = "";
//         for( unsigned int rng_idx = 0;
//             ( rng_idx < rangesFound.size() );
//             ++rng_idx )
//         {
//             out << pos_prefix;
//             const TRange &range = rangesFound[rng_idx];
//             out << (range.first + 1); // "+1" because 1-based for user
// ```
// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta_exception.hpp:90-92
// ```c++
//         void ConvertBadIndexesToString(
//             CNcbiOstream & out,
//             unsigned int maxRanges = 1000 ) const;
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2117-2125
// ```c++
// std::string CFastaReader::x_NucOrProt(void) const
// {
//     if( m_CurrentSeq && m_CurrentSeq->IsSetInst() &&
//         m_CurrentSeq->GetInst().IsSetMol() )
//     {
//         return ( m_CurrentSeq->GetInst().IsAa() ? "protein " : "nucleotide " );
//     } else {
//         return kEmptyStr;
//     }
// ```
// NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:226-230
// ```c++
// #define LOG_POST(message)                                               \
//     ( NCBI_NS_NCBI::CNcbiDiag(DIAG_COMPILE_INFO,                        \
//       NCBI_NS_NCBI::eDiag_Error,                                        \
//       NCBI_NS_NCBI::eDPF_Log | NCBI_NS_NCBI::eDPF_IsNote).GetRef()      \
//       << message                                                        \
// ```
// NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:747-747
// ```c++
//     eDPF_Log                = 0,
// ```
// ParseDataLine precedes AssignMolType (fasta.cpp:409,1360), so fresh
// BLAST FASTA records have no molecule label when this warning is posted.
// LOG_POST_X uses the log flags without severity/prefix. Positions in the
// existing B reader state have already been converted to one-based indexes.
pub(crate) fn reader_warning(warning: &InputWarning, stderr: &mut dyn Write) -> Result<()> {
    let message=match warning.kind {
            "hyphens"=>format!("CFastaReader: Hyphens are invalid and will be ignored around line {}",warning.line),
            "title_nucleotide_residues"=>"FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?".into(),
            "title_amino_acids"=>"FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?".into(),
            "invalid_residues"=>{
                let mut ranges:Vec<(usize,usize)>=Vec::new();
                for &p in &warning.positions {if let Some(last)=ranges.last_mut().filter(|r|r.1+1==p){last.1=p;}else{if ranges.len()>=1000{break;}ranges.push((p,p));}}
                let text=ranges.iter().map(|&(a,b)|if a==b{a.to_string()}else{format!("{a}-{b}")}).collect::<Vec<_>>().join(", ");
                format!("FASTA-Reader: Ignoring invalid {}residues at position(s): On line {}: {text}",String::from_utf8_lossy(&warning.text),warning.line)
            },
            "defline_id"=>String::from_utf8_lossy(&warning.text).into_owned(),
            _=>unreachable!("source-registered FASTA warning"),
        };
    writeln!(stderr, "{message}")?;
    Ok(())
}

/// NCBI searches a reader batch as a whole: it is not searched only when none of its contexts is
/// valid. The units of a batch are searched on their own, so a unit without a valid context is
/// marked skipped; `search_skipped` of all of them is set from the whole batch.
fn reader_batch_skipped(results: &mut [BatchResults]) {
    let skipped = results.iter().all(|results| results.search_skipped);
    for results in results {
        results.search_skipped = skipped;
    }
}

// No NCBI counterpart: how LOSAT schedules the search of a reader batch (NCBI's query batch,
// `BATCH_SIZE` letters). The search keeps units of at most 10002 letters (`batch_ranges` with
// the constant), whatever `BATCH_SIZE` is: the flow and the waves schedule the units as they
// scheduled the batches of a run without `BATCH_SIZE`, and the report of a reader batch is
// written when its last unit is searched. The values that NCBI takes of the whole query batch
// (the split decision, the ungapped link cutoffs, the sum-statistics average length, the
// composition redo's minimum cutoff) are taken of the reader batch, not of the unit: phase 1
// prepares every unit of the batch and reduces them (`reader_batch.rs`), phase 2 searches the
// units with the result. The other values NCBI takes of the batch (the first valid context, the
// search-skipped flag) are taken by `ReaderBatchResults` below.
/// One unit of the search: `range` of a reader batch, which `last` ends.
#[cfg(feature = "parallel")]
#[derive(Clone)]
struct XUnit {
    batch: std::sync::Arc<Vec<FastaRecord>>,
    range: std::ops::Range<usize>,
    last: bool,
    /// Phase 1 of the unit's reader batch when the batch has several units (or, in the capped batch
    /// flow, `LOSAT_X_BXCAPFILL`, is the unit itself), and the unit's place in it. `None`: the unit
    /// is the whole reader batch and prepares itself when it is searched.
    prep: Option<(std::sync::Arc<XReaderPrep>, usize)>,
}

/// Phase 1 of a reader batch of several units (`reader_batch.rs`), or of one unit prepared ahead
/// (`LOSAT_X_BXCAPFILL`): the units' preparations, and (the flow) the searches that wait for them.
#[cfg(feature = "parallel")]
struct XReaderPrep {
    /// The batch's values that follow from its records; each unit is prepared with them. `None`:
    /// a reader batch of one unit prepared ahead of its search (`LOSAT_X_BXCAPFILL`), which has no
    /// shape (`reader_batch::prepare_lone_unit`).
    layout: Option<std::sync::Arc<ReaderBatchShape>>,
    /// The units' preparations, from the end of phase 1 until each unit's search takes its own.
    /// All `None` if phase 1 failed: each unit then prepares itself when it is searched, and the
    /// one whose queries fail raises the error there (a reader batch with a search error writes
    /// nothing, so the values of the units before it are not used).
    slots: Vec<std::sync::Mutex<Option<UnitPrep>>>,
    /// The flow's tasks of phase 1 not yet ended, and whether one failed.
    remaining: std::sync::atomic::AtomicUsize,
    failed: std::sync::atomic::AtomicBool,
    /// Whether phase 1 has ended, and the units whose search was due before it ended (the flow:
    /// number, unit, searched alone).
    gate: std::sync::Mutex<(bool, Vec<(usize, XUnit, bool)>)>,
}

#[cfg(feature = "parallel")]
impl XReaderPrep {
    /// The end of phase 1: the batch's shape is reduced and attached to every preparation (or,
    /// if a unit failed, no unit keeps one), and the searches that waited are handed back.
    fn finish(&self) -> Vec<(usize, XUnit, bool)> {
        use std::sync::atomic::Ordering;
        let mut slots: Vec<_> = self
            .slots
            .iter()
            .map(|slot| slot.lock().expect("a unit's preparation"))
            .collect();
        if self.failed.load(Ordering::Acquire) {
            slots.iter_mut().for_each(|slot| **slot = None);
        } else if let Some(layout) = &self.layout {
            let shape = std::sync::Arc::new(ReaderBatchShape::reduce(
                layout,
                slots
                    .iter()
                    .map(|slot| slot.as_ref().expect("a prepared unit")),
            ));
            for slot in &mut slots {
                attach(slot.as_mut().expect("a prepared unit"), &shape);
            }
        }
        drop(slots);
        let mut gate = self.gate.lock().expect("the reader batch's gate");
        gate.0 = true;
        std::mem::take(&mut gate.1)
    }
}

#[cfg(feature = "parallel")]
impl XUnit {
    /// The units of a reader batch, in order (the last one marked); a copy of the batch is
    /// shared by them. A batch of several units gets its phase 1 (`prepare_now`, or the flow's
    /// tasks), which its units' searches take; with `ahead` (the capped batch flow,
    /// `LOSAT_X_BXCAPFILL`) a batch of one unit of at most two search units' letters gets one too
    /// (the flow's task, without a shape). A longer batch (one long query) prepares itself when it
    /// is searched, as before: preparations held ahead are bounded by letters, not only by the
    /// flow's count of unwritten batches.
    fn split(
        queries: &[FastaRecord],
        options: &super::args::ResolvedOptions,
        ahead: bool,
    ) -> Vec<XUnit> {
        let batch = std::sync::Arc::new(queries.to_vec());
        let ranges = super::query_setup::batch_ranges(&batch, BATCH_SIZE);
        let count = ranges.len();
        let ahead = ahead
            && count == 1
            && batch.iter().map(|q| q.sequence.len()).sum::<usize>() <= 2 * BATCH_SIZE;
        let prep = (count > 1 || ahead).then(|| {
            std::sync::Arc::new(XReaderPrep {
                layout: (count > 1)
                    .then(|| std::sync::Arc::new(ReaderBatchShape::layout(&batch, options))),
                slots: (0..count).map(|_| std::sync::Mutex::new(None)).collect(),
                remaining: std::sync::atomic::AtomicUsize::new(count),
                failed: std::sync::atomic::AtomicBool::new(false),
                gate: std::sync::Mutex::new((false, Vec::new())),
            })
        });
        ranges
            .into_iter()
            .enumerate()
            .map(|(index, range)| XUnit {
                batch: std::sync::Arc::clone(&batch),
                range,
                last: index + 1 == count,
                prep: prep
                    .as_ref()
                    .map(|prep| (std::sync::Arc::clone(prep), index)),
            })
            .collect()
    }

    /// Phase 1 of the units' reader batch (`units` from `split`), all at once by `prepare`.
    fn prepare_now(
        units: &[XUnit],
        prepare: &dyn Fn(&[FastaRecord], &[std::ops::Range<usize>]) -> Result<Vec<UnitPrep>>,
    ) {
        let Some((prep, _)) = units.first().and_then(|unit| unit.prep.as_ref()) else {
            return;
        };
        let ranges: Vec<_> = units.iter().map(|unit| unit.range.clone()).collect();
        if let Ok(preps) = prepare(&units[0].batch, &ranges) {
            for (slot, unit) in prep.slots.iter().zip(preps) {
                *slot.lock().expect("a unit's preparation") = Some(unit);
            }
        }
        prep.gate.lock().expect("the reader batch's gate").0 = true;
    }

    /// The flow's task of phase 1 for this unit: prepares it by `prepare` (with the batch's layout,
    /// or none for a batch of one unit prepared ahead); the last task of the reader batch ends
    /// phase 1 and returns the searches that waited for it.
    fn prepare_task(
        &self,
        prepare: &(dyn Fn(&[FastaRecord], Option<&std::sync::Arc<ReaderBatchShape>>) -> Result<UnitPrep>
              + Sync),
    ) -> Vec<(usize, XUnit, bool)> {
        use std::sync::atomic::Ordering;
        let Some((prep, position)) = self.prep.as_ref() else {
            return Vec::new();
        };
        // A panic is a failure of phase 1: the unit prepares itself again when it is searched,
        // and its search raises the panic where the unit is written. (With `panic = "abort"`, the
        // release profile, the process ends at the first panic instead, which may be before the
        // batches ahead of this one are written.)
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            prepare(self.queries(), prep.layout.as_ref())
        })) {
            Ok(Ok(unit)) => {
                *prep.slots[*position].lock().expect("a unit's preparation") = Some(unit)
            }
            _ => prep.failed.store(true, Ordering::Release),
        }
        if prep.remaining.fetch_sub(1, Ordering::AcqRel) == 1 {
            prep.finish()
        } else {
            Vec::new()
        }
    }

    /// The unit, if its search can start now: a unit of a reader batch whose phase 1 has not
    /// ended waits for it (number `index`, searched `alone`), and its search is started when
    /// phase 1 ends.
    fn ready(self, index: usize, alone: bool) -> Option<XUnit> {
        let Some((prep, _)) = self.prep.as_ref() else {
            return Some(self);
        };
        let prep = std::sync::Arc::clone(prep);
        let mut gate = prep.gate.lock().expect("the reader batch's gate");
        if gate.0 {
            return Some(self);
        }
        gate.1.push((index, self, alone));
        None
    }

    /// The unit's phase 1 preparation, once.
    fn take_prep(&self) -> Option<UnitPrep> {
        self.prep.as_ref().and_then(|(prep, position)| {
            prep.slots[*position]
                .lock()
                .expect("a unit's preparation")
                .take()
        })
    }

    fn queries(&self) -> &[FastaRecord] {
        &self.batch[self.range.clone()]
    }
}

/// The results of the units of a reader batch searched so far. At the batch's last unit they are
/// handed over as the results of the reader batch: `query_ordinal` counts from the reader
/// batch's first query, and `search_skipped` is set for all of them when no unit had a valid
/// context (NCBI searches the batch, or does not, as a whole).
#[cfg(feature = "parallel")]
#[derive(Default)]
struct ReaderBatchResults {
    results: Vec<BatchResults>,
}

#[cfg(feature = "parallel")]
impl ReaderBatchResults {
    fn add(
        &mut self,
        start: usize,
        last: bool,
        mut unit: Vec<BatchResults>,
    ) -> Option<Vec<BatchResults>> {
        for results in &mut unit {
            results.query_ordinal += start;
            // The report reads the contexts of the prepared queries and nothing of their
            // residues, which are not needed once the unit is searched; the results of a unit
            // wait here for the last unit of the batch.
            results.prepared.sequence_start = Vec::new();
            results.prepared.sequence_start_nomask = Vec::new();
            results.prepared.lookup_segments = Vec::new();
        }
        self.results.append(&mut unit);
        if !last {
            return None;
        }
        let mut results = std::mem::take(&mut self.results);
        reader_batch_skipped(&mut results);
        Some(results)
    }
}

/// A unit's search in the flow, or the panic it ended with (raised again where it is written).
#[cfg(feature = "parallel")]
type XSearched = std::thread::Result<Result<Vec<BatchResults>>>;

// EXPERIMENT (LOSAT_X_BXFLOW): the reading loop of `run` with every query batch searched as a task
// of its own (see the call in `run`). `search` searches one batch (the waves' function; `alone`:
// with the pool inside), `write_batch` writes one batch's reports (or ends the run at its error),
// `write_warning` writes the bytes of one reader warning. A batch is spawned once `lookahead`
// batches follow it and while fewer than `max_live` are spawned and not yet searched (by the
// reading thread or by a search that ends); the batches still held when the input ends are
// searched alone. At most `bound` batches are read and not yet written. `split` makes a reader
// batch its units (the flow's batches); a reader batch of several units (in the capped flow with
// `LOSAT_X_BXCAPFILL`, every reader batch) is prepared by `prepare` (phase 1, one task per unit)
// before any of its units is searched. Returns what the reader returns.
// No NCBI counterpart: scheduling of the batch searches only (NCBI's -mt_mode 1 runs batches as
// nodes and prints them in order, blastx_app.cpp:345-357, cited in `run`); every value written
// is computed by the same calls as without the switch, in the same order of writing.
#[cfg(feature = "parallel")]
#[allow(clippy::too_many_arguments)]
fn x_flow(
    bound: usize,
    lookahead: usize,
    max_live: usize,
    query_stream: &mut FastaStream<File>,
    lcase_masking: bool,
    batch_size: usize,
    keep_defline: bool,
    split: &dyn Fn(&[FastaRecord]) -> Vec<XUnit>,
    prepare: &(dyn Fn(&[FastaRecord], Option<&std::sync::Arc<ReaderBatchShape>>) -> Result<UnitPrep>
          + Sync),
    search: &(dyn Fn(&XUnit, bool) -> Result<Vec<BatchResults>> + Sync),
    write_unit: &mut dyn FnMut(
        &[FastaRecord],
        usize,
        bool,
        Result<Vec<BatchResults>>,
    ) -> Result<()>,
    write_warning: &mut dyn FnMut(&[u8]) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    use super::x_timing::Run as XRun;
    use std::cell::{Cell, RefCell};
    use std::collections::{BTreeMap, VecDeque};
    use std::sync::mpsc;
    // A reader event, written in the reader's order: a warning's bytes, or a unit by its number.
    enum XEvent {
        Warn(Vec<u8>),
        Batch(usize, XUnit),
    }
    struct XState {
        events: VecDeque<XEvent>,
        searched: BTreeMap<usize, XSearched>,
        read: usize,
        unwritten: usize,
    }
    // Units read and not yet spawned (the newest `lookahead` until the end of the input), in
    // reader order, and the number spawned and not yet searched (at most `max_live`; a unit that
    // waits for its reader batch's phase 1 is spawned).
    struct XHeld {
        batches: VecDeque<(usize, XUnit)>,
        live: usize,
        ended: bool,
        // The writing ended (at the end, at an error, or at a panic it raises again): spawn no more.
        stopped: bool,
    }
    struct XLaunch<'f> {
        held: std::sync::Mutex<XHeld>,
        lookahead: usize,
        max_live: usize,
        search: &'f (dyn Fn(&XUnit, bool) -> Result<Vec<BatchResults>> + Sync),
    }
    impl XLaunch<'_> {
        // The held units; a panic elsewhere cannot leave them unusable (nothing panics while
        // they are locked, and a search's result must still be sent).
        fn held(&self) -> std::sync::MutexGuard<'_, XHeld> {
            self.held
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
        }

        // Spawns the held units beyond the newest `lookahead` (all of them at the end of the
        // input, then searched alone) while fewer than `max_live` are spawned and unsearched. A
        // unit of a reader batch whose phase 1 has not ended is searched when it ends.
        fn release<'s>(
            &'s self,
            scope: &rayon::Scope<'s>,
            sender: &mpsc::Sender<(usize, XSearched)>,
        ) {
            loop {
                let next = {
                    let mut held = self.held();
                    let keep = if held.ended { 0 } else { self.lookahead };
                    if !held.stopped && held.batches.len() > keep && held.live < self.max_live {
                        held.live += 1;
                        let alone = held.ended;
                        held.batches
                            .pop_front()
                            .map(|(index, unit)| (index, unit, alone))
                    } else {
                        None
                    }
                };
                let Some((index, unit, alone)) = next else {
                    break;
                };
                if alone {
                    super::x_timing::run_add(XRun::FlowAlone, 1);
                }
                if let Some(unit) = unit.ready(index, alone) {
                    self.spawn_search(scope, sender, index, unit, alone);
                }
            }
        }

        // Spawns the search of unit `index`, spawned by `release`. The search sends its result
        // (or panic) once and then spawns again, so a place freed by a search is filled without
        // waiting for the writing thread. A unit whose phase 1 ends after the writing ended is
        // not searched.
        fn spawn_search<'s>(
            &'s self,
            scope: &rayon::Scope<'s>,
            sender: &mpsc::Sender<(usize, XSearched)>,
            index: usize,
            unit: XUnit,
            alone: bool,
        ) {
            {
                let mut held = self.held();
                if held.stopped {
                    held.live -= 1;
                    return;
                }
            }
            let sender = sender.clone();
            scope.spawn(move |scope| {
                let searched = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    (self.search)(&unit, alone)
                }));
                self.held().live -= 1;
                // The receiver is gone only when the run already ended at an error.
                let _ = sender.send((index, searched));
                self.release(scope, &sender);
            });
        }
    }
    debug_assert!(rayon::current_thread_index().is_some());
    debug_assert!(lookahead < bound);
    let launch = XLaunch {
        held: std::sync::Mutex::new(XHeld {
            batches: VecDeque::new(),
            live: 0,
            ended: false,
            stopped: false,
        }),
        lookahead,
        max_live,
        search,
    };
    // Marks the writing ended when the scope's body returns or unwinds, before the scope waits
    // for the searches it spawned: the searches that end then spawn no more.
    struct XStop<'a, 'f>(&'a XLaunch<'f>);
    impl Drop for XStop<'_, '_> {
        fn drop(&mut self) {
            self.0.held().stopped = true;
        }
    }
    rayon::in_place_scope(|scope| {
        let _stop = XStop(&launch);
        let (sender, receiver) = mpsc::channel::<(usize, XSearched)>();
        let state = RefCell::new(XState {
            events: VecDeque::new(),
            searched: BTreeMap::new(),
            read: 0,
            unwritten: 0,
        });
        let failed = Cell::new(false);
        let release = || launch.release(scope, &sender);
        // Writes the events at the front whose batches are searched; then, until fewer than
        // `bound` batches are unwritten (or, with `all`, until every event is written), runs a
        // pending task or waits for a result, and writes again.
        let mut write = |all: bool| -> Result<()> {
            loop {
                {
                    let mut state = state.borrow_mut();
                    while let Ok((index, searched)) = receiver.try_recv() {
                        state.searched.insert(index, searched);
                    }
                }
                loop {
                    let event = {
                        let mut state = state.borrow_mut();
                        let ready = match state.events.front() {
                            Some(XEvent::Warn(_)) => true,
                            Some(XEvent::Batch(index, _)) => state.searched.contains_key(index),
                            None => false,
                        };
                        if ready {
                            state.events.pop_front()
                        } else {
                            None
                        }
                    };
                    let x_t = super::x_timing::lap();
                    match event {
                        None => break,
                        Some(XEvent::Warn(bytes)) => write_warning(&bytes)?,
                        Some(XEvent::Batch(index, unit)) => {
                            let searched = {
                                let mut state = state.borrow_mut();
                                state.unwritten -= 1;
                                state.searched.remove(&index).expect("a searched batch")
                            };
                            match searched {
                                Ok(results) => {
                                    write_unit(&unit.batch, unit.range.start, unit.last, results)?
                                }
                                Err(panic) => std::panic::resume_unwind(panic),
                            }
                        }
                    }
                    super::x_timing::run(XRun::Render, x_t);
                }
                {
                    let state = state.borrow();
                    if (all && state.events.is_empty()) || (!all && state.unwritten < bound) {
                        return Ok(());
                    }
                }
                let x_t = super::x_timing::lap();
                if let Some(rayon::Yield::Executed) = rayon::yield_now() {
                    let x_ns = super::x_timing::run(XRun::FlowJobTime, x_t);
                    super::x_timing::run_add(XRun::FlowJobs, 1);
                    super::x_timing::run_add(XRun::Search, x_ns);
                } else {
                    // No task is pending now: the unwritten batches are searched on other
                    // threads (each sends its result once). Wait for a result, but look for
                    // pending work again soon: a batch searched alone hands out parts of its
                    // search as it goes.
                    match receiver.recv_timeout(std::time::Duration::from_millis(1)) {
                        Ok((index, searched)) => {
                            state.borrow_mut().searched.insert(index, searched);
                        }
                        Err(mpsc::RecvTimeoutError::Timeout) => {}
                        Err(mpsc::RecvTimeoutError::Disconnected) => {
                            unreachable!("the reading thread keeps a sender")
                        }
                    }
                    let x_ns = super::x_timing::run(XRun::FlowBlocked, x_t);
                    super::x_timing::run_add(XRun::Search, x_ns);
                }
            }
        };
        let parsed = parse_fasta_stream_keeping(
            query_stream,
            false,
            lcase_masking,
            batch_size,
            keep_defline,
            &mut |queries| {
                if queries.is_empty() {
                    return Err(engine("Empty CBlastQueryVector"));
                }
                // The reader batch enters the flow as its units (at most 10002 letters each),
                // each as a batch of a run without `BATCH_SIZE` does; the bound counts units. A
                // reader batch of several units (in the capped flow with LOSAT_X_BXCAPFILL, every
                // reader batch) is prepared first (phase 1) by one task per unit,
                // beside the searches of the batches before it, while this thread reads on: the
                // last task ends phase 1 and starts the searches that were due before it ended
                // (no thread waits for phase 1).
                let units = split(queries);
                if units.first().is_some_and(|unit| unit.prep.is_some()) {
                    for unit in &units {
                        let (unit, sender, launch) = (unit.clone(), sender.clone(), &launch);
                        scope.spawn(move |scope| {
                            for (index, unit, alone) in unit.prepare_task(prepare) {
                                launch.spawn_search(scope, &sender, index, unit, alone);
                            }
                        });
                    }
                }
                for unit in units {
                    {
                        let mut state = state.borrow_mut();
                        let index = state.read;
                        state.read += 1;
                        state.unwritten += 1;
                        state.events.push_back(XEvent::Batch(index, unit.clone()));
                        let mut held = launch.held();
                        held.batches.push_back((index, unit));
                    }
                    super::x_timing::run_add(XRun::FlowBatches, 1);
                    release();
                    let written = write(false);
                    if written.is_err() {
                        failed.set(true);
                        return written;
                    }
                }
                Ok(())
            },
            &mut |warning| {
                let mut bytes = Vec::new();
                reader_warning(warning, &mut bytes)?;
                state.borrow_mut().events.push_back(XEvent::Warn(bytes));
                Ok(())
            },
        );
        if !failed.get() {
            // The end of the input: the held batches (the last of the input) are searched alone,
            // side by side.
            launch.held().ended = true;
            release();
            write(true)?;
        }
        parsed
    })
}

/// The formatting handler's warnings of a BLASTX run with output format `outfmt`
/// (`ReportCounts::warnings`).
pub(crate) fn warnings_of(args: &BlastxArgs, outfmt: i32) -> Vec<u8> {
    crate::blastinput::report_counts::report_counts(
        outfmt,
        args.max_target_seqs.map(|n| n as usize),
        args.num_descriptions,
        args.num_alignments,
    )
    .warnings("blastx", outfmt, args.num_descriptions.is_some())
}

#[cfg(test)]
mod archive_error_tests {
    use super::*;
    use crate::report::serial::encoder::SerialError;

    /// The encoder's refusal stays LOSAT's own error: no `BLAST engine error:` text and no
    /// error archive. NCBI's engine errors still become both.
    #[test]
    fn the_encoders_refusal_gets_no_error_archive() {
        let refusal = archive_error(SerialError::NonVisibleChar(0xc3).into());
        assert!(refusal.is::<SerialError>());
        assert!(archive::error_lines(&refusal, false).is_none());
        assert!(archive::error_lines(&refusal, true).is_none());
        let engine_error = archive_error(anyhow::anyhow!("something else"));
        assert!(engine_error.downcast_ref::<NativeError>().is_some());
        assert!(archive::error_lines(&engine_error, false).is_some());
    }
}

#[cfg(all(test, feature = "parallel"))]
mod reader_batch_tests {
    use super::*;

    fn record(letters: usize) -> FastaRecord {
        FastaRecord {
            internal_id: String::new(),
            title: String::new(),
            sequence: vec![b'A'; letters],
            lowercase_masks: Vec::new(),
            warnings: Vec::new(),
            source: None,
        }
    }

    /// One unit's result for `queries` queries, skipped or not.
    fn unit_results(queries: usize, skipped: bool) -> Vec<BatchResults> {
        let records: Vec<_> = (0..queries).map(|_| record(30)).collect();
        vec![BatchResults {
            query_ordinal: 0,
            queries: vec![Vec::new(); queries],
            prepared: super::super::query_setup::prepare_queries(
                &records,
                &super::super::stage_c_tests::options(),
            )
            .unwrap(),
            parameters: Vec::new(),
            search_skipped: skipped,
        }]
    }

    /// A reader batch (any `BATCH_SIZE`) is split into units of at most 10002 letters, as a run
    /// without `BATCH_SIZE` batches its queries; the last unit is marked, and the units cover the
    /// batch in order.
    #[test]
    fn a_reader_batch_is_split_into_units_of_10002_letters() {
        let batch = [4000, 4000, 4000, 100, 10002, 5, 20000, 1]
            .map(record)
            .to_vec();
        let options = super::super::stage_c_tests::options();
        let units = XUnit::split(&batch, &options, false);
        let ranges: Vec<_> = units.iter().map(|unit| unit.range.clone()).collect();
        // 4000 + 4000 + 4000 reaches 10002 at the third record; 100 + 10002 at the fifth;
        // 5 + 20000 at the seventh; the rest is the last unit.
        assert_eq!(ranges, vec![0..3, 3..5, 5..7, 7..8]);
        assert_eq!(
            units.iter().map(|unit| unit.last).collect::<Vec<_>>(),
            [false, false, false, true]
        );
        assert_eq!(units[1].queries().len(), 2);
        // A reader batch of at most 10002 letters (the default, or a smaller `BATCH_SIZE`) is
        // one unit.
        let small = XUnit::split(&batch[..3], &options, false);
        assert_eq!(small.len(), 1);
        assert!(small[0].last);
    }

    /// The results of the units become the reader batch's at its last unit: ordinals count from
    /// the batch's first query, and a unit with no valid context is skipped only when no unit
    /// had one.
    #[test]
    fn the_units_results_are_the_reader_batchs() {
        let mut pending = ReaderBatchResults::default();
        assert!(pending.add(0, false, unit_results(3, false)).is_none());
        assert!(pending.add(3, false, unit_results(1, true)).is_none());
        let results = pending.add(4, true, unit_results(2, false)).unwrap();
        let ordinals: Vec<_> = results.iter().map(|r| r.query_ordinal).collect();
        assert_eq!(ordinals, [0, 3, 4]);
        assert!(results.iter().all(|r| !r.search_skipped));
        // Nothing is left for the next reader batch, and a batch of skipped units is skipped.
        let results = pending
            .add(0, false, unit_results(1, true))
            .or_else(|| pending.add(1, true, unit_results(1, true)))
            .unwrap();
        assert_eq!(results.len(), 2);
        assert!(results.iter().all(|r| r.search_skipped));
    }
}
