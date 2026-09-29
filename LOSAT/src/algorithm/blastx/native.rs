//! Native serial caller, tested through the internal probe before public dispatch.
use super::{
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
        parse_fasta_stream_with_warnings, stream_is_empty, FastaParseError, FastaStream,
        InputWarning,
    },
    query_setup::BATCH_SIZE,
    report,
    runtime::search_internal,
};
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
pub use crate::cli::{exit_on_native_error, NativeError};
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:95-99
// ```c++
// string s_ArgExptMsg(const string& name, const string& what, const string& attr)
// {
//     return string("Argument \"") + (name.empty() ? s_ExtraName : name) +
//         "\". " + what + (attr.empty() ? attr : ":  `" + attr + "'");
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-619
// ```c++
// void CArg_Ios::x_Open(CArgValue::TFileFlags /*flags*/) const
// {
//     if ( !m_Ios ) {
//         NCBI_THROW(CArgException,eNoFile, s_ArgExptMsg(GetName(),
//             "File is not accessible",AsString()));
// ```
fn inaccessible(name: &str, path: &Path) -> anyhow::Error {
    NativeError {
        exit: 1,
        message: format!(
            "Command line argument error: Argument \"{name}\". File is not accessible:  `{path}'\n",
            path = path.display()
        ),
    }
    .into()
}
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
    if error.downcast_ref::<report::TabularFlushError>().is_some() {
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
pub fn run(args: &BlastxArgs) -> Result<()> {
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
    let subject_file =
        File::open(&args.subject).map_err(|_| inaccessible("subject", &args.subject))?;
    let mut subject_stream = FastaStream::from_file(subject_file);
    let mut stderr = io::stderr().lock();
    let subjects = parse_fasta_stream_with_warnings(
        &mut subject_stream,
        true,
        args.lcase_masking,
        usize::MAX,
        &mut |_| Ok(()),
        &mut |warning| reader_warning(warning, &mut stderr),
    )
    .map_err(input_error)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    // ```c++
    // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
    //     : m_QueryVector(& queries)
    // {
    //     if (queries.Empty()) {
    //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
    //     }
    // ```
    if subjects.is_empty() {
        return Err(engine("Empty CBlastQueryVector"));
    }
    let query_file = File::open(&args.query).map_err(|_| inaccessible("query", &args.query))?;
    let destination: Box<dyn Write> = match args.out.as_deref().filter(|p| *p != Path::new("-")) {
        Some(path) => Box::new(File::create(path).map_err(|_| inaccessible("out", path))?),
        None => Box::new(io::stdout().lock()),
    };
    let mut writer = BufWriter::new(destination);
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
    // ```c++
    //     if(hitlist_size < 5){
    //    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
    //     }
    // ```
    if args.max_target_seqs.unwrap_or(500) < 5 {
        stderr.write_all(b"Warning: [blastx] Examining 5 or more matches is recommended\n")?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-175
    // ```c++
    //     catch (const blast::CInputException& e) {                               \
    //         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
    //         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
    //         exit_code = BLAST_INPUT_ERROR;                                      \
    // ```
    let options = resolution.map_err(|e| NativeError {
        exit: 1,
        message: format!(
            "BLAST query/options error: {e}\nPlease refer to the BLAST+ user manual.\n"
        ),
    })?;
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
    let mut query_stream = FastaStream::from_file(query_file);
    if stream_is_empty(&mut query_stream) {
        stderr.write_all(b"Warning: [blastx] Query is Empty!\n")?;
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:217-218
    // ```c++
    //         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
    //         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
    // ```
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
            s.internal_id, s.title
        )?;
    }
    report::write_prolog(&mut writer, &options, &subjects, &args.subject).map_err(report_error)?;
    let mut processed = 0usize;
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
        parse_fasta_stream_with_warnings(
            &mut query_stream,
            false,
            args.lcase_masking,
            BATCH_SIZE,
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
                let mut results = search_internal(queries, &subjects, &options).map_err(engine)?;
                report::render_queries(
                    &mut writer,
                    &mut **diagnostics.borrow_mut(),
                    &options,
                    queries,
                    &subjects,
                    &args.subject,
                    &mut results,
                )
                .map_err(report_error)?;
                writer.flush().map_err(output_error)?;
                processed += queries.len();
                Ok(())
            },
            &mut |warning| reader_warning(warning, &mut **diagnostics.borrow_mut()),
        )
        .map_err(input_error)?;
    }
    report::write_epilog(&mut writer, &options, &subjects, &args.subject, processed)
        .map_err(report_error)?;
    writer.flush().map_err(output_error)?;
    Ok(())
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
                format!("FASTA-Reader: Ignoring invalid residues at position(s): On line {}: {text}",warning.line)
            },_=>unreachable!("source-registered FASTA warning"),
        };
    writeln!(stderr, "{message}")?;
    Ok(())
}
