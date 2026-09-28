//! Local FASTA state before translation. Internal IDs and defline titles differ.
use anyhow::{bail, Context, Result};
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
// ```c++
// bool CStreamLineReader::AtEOF(void) const
// {
//     return !m_UngetLine &&
//         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
// }
// ```
use std::{
    io::{BufRead, BufReader, Read, Seek, SeekFrom},
    path::Path,
};

#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:619-631
// ```c++
// void CFastaReader::PostProcessIDs(
//     const CBioseq::TId& defline_ids,
//     const string& /*defline*/,
//     const bool has_range,
//     const TSeqPos range_start,
//     const TSeqPos range_end)
// {
//     if (defline_ids.empty()) {
//         GenerateID();
//     }
//     else {
//         SetIDs() = defline_ids;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:177-180
// ```c++
//     size_t title_start = NPOS;
//     if ((fFastaFlags & CFastaReader::fNoParseID)) {
//         title_start = start;
//     }
// ```
pub struct FastaRecord {
    pub internal_id: String,
    pub title: String,
    pub sequence: Vec<u8>,
    pub lowercase_masks: Vec<(i32, i32)>,
    pub warnings: Vec<InputWarning>,
}
#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:974-995
// ```c++
//         case eCharType_Bad:
//             if( bad_pos_line_num < 0 ) {
//                 bad_pos_line_num = LineNumber();
//             }
//             bad_pos_vec.push_back(pos);
//             break;
//         default:
//             _TROUBLE;
//         }
//     }
//
//     m_SeqData.resize(m_CurrentPos);
//
//     if( bIgnorableHyphenSeen ) {
//         _ASSERT( bHyphensIgnoreAndWarn );
//         FASTA_WARNING_EX(LineNumber(),
//             "CFastaReader: Hyphens are invalid and will be ignored around line " << LineNumber(),
//             ILineError::eProblem_IgnoredResidue,
//             kEmptyStr, kEmptyStr, "-" );
//     }
//
//     // before throwing, be sure that we're in a valid state so that callers can
// ```
pub struct InputWarning {
    pub line: usize,
    pub kind: &'static str,
    pub positions: Vec<usize>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:319-343
// ```c++
//     CFastaReader::TFlags flags = m_Config.GetBelieveDeflines() ?
//                                     CFastaReader::fParseRawID:
//                                     (CFastaReader::fNoParseID |
//                                      CFastaReader::fDLOptional);
//
//     // Allow CFastaReader fSkipCheck flag to be set based
//     // on new CBlastInputSourceConfig property - GetSkipSeqCheck() -RMH-
//     flags += ( m_Config.GetSkipSeqCheck() ? CFastaReader::fSkipCheck : 0 );
//
//     flags += (m_ReadProteins
//               ? CFastaReader::fAssumeProt
//               : CFastaReader::fAssumeNuc);
//     const char* env_var = getenv("BLASTINPUT_GEN_DELTA_SEQ");
//     if (env_var == NULL || (env_var && string(env_var) == kEmptyStr)) {
//         flags += CFastaReader::fNoSplit;
//     }
//     // This is necessary to enable the ignoring of gaps in classes derived from
//     // CFastaReader
//
//    	flags+= CFastaReader::fHyphensIgnoreAndWarn;
//
//     flags+= CFastaReader::fDisableNoResidues;
//     // Do not check more than few characters in local ID for illegal characters.
//     // Illegal characters can be things like = and we want to let those through.
//     flags+= CFastaReader::fQuickIDCheck;
// ```
pub fn read_fasta(path: &Path, protein: bool, lowercase: bool) -> Result<Vec<FastaRecord>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:291-301
    // ```c++
    // CBlastFastaInputSource::CBlastFastaInputSource(CNcbiIstream& infile,
    //                                        const CBlastInputSourceConfig& iconfig)
    //     : m_Config(iconfig),
    //       m_LineReader(iconfig.GetConvertGapsToNs() ?
    //                    new CStreamLineReaderConverter(infile) :
    //                    new CStreamLineReader(infile)),
    //       m_ReadProteins(iconfig.IsProteinInput())
    // {
    //     x_InitInputReader();
    // }
    //
    // ```
    let file = std::fs::File::open(path)
        .with_context(|| format!("cannot read BLASTX input {}", path.display()))?;
    parse_fasta_stream_with_warnings(
        &mut FastaStream::from_file(file),
        protein,
        lowercase,
        usize::MAX,
        &mut |_| Ok(()),
        &mut |_| Ok(()),
    )
}
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:350-375
// ```c++
//         if (c == '>' ) {
//             CTempString next_line = *++GetLineReader();
//             string strmodified;
//             if( NStr::StartsWith(next_line, ">?_") ) {
//                 CTempString tmp = next_line.substr(3);
//                 strmodified = ">";
//                 strmodified.append(tmp.data(), tmp.length());
//                 next_line = strmodified;
//             }
//             if( NStr::StartsWith(next_line, ">?") ) {
//                 // This is actually a data line. an assembly gap, in particular, which
//                 // we handle farther below
//                 GetLineReader().UngetLine();
//             } else {
//                 if(need_defline) {
//                     ParseDefLine(next_line, pMessageListener);
//                     need_defline = false;
//                     continue;
//                 } else {
//                     GetLineReader().UngetLine();
//                     // start of the next sequence
//                     break;
//                 }
//             }
//         }
//
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:378-389
// ```c++
//         if (line.empty()) {
//             continue; // ignore lines containing only whitespace
//         }
//         c = line[0];
//
//         if (c == '!'  ||  c == '#' || c == ';') {
//             // no content, just a comment or blank line
//             continue;
//         } else if (need_defline) {
//             if (TestFlag(fDLOptional)) {
//                 ParseDefLine(">", pMessageListener);
//                 need_defline = false;
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:941-966
// ```c++
//             m_SeqData[m_CurrentPos] = c;
//             CloseMask();
//             ++m_CurrentPos;
//             break;
//         case eCharType_MaskedNonGap:
//             CloseGap(pos == 0);
//             m_SeqData[m_CurrentPos] = s_ASCII_MustBeLowerToUpper(c);
//             OpenMask();
//             ++m_CurrentPos;
//             break;
//         case eCharType_Gap: {
//             CloseMask();
//             // open a gap
//
//             size_t pos2 = pos + 1;
//             while( pos2 < s_len && s[pos2] == c ) {
//                 ++pos2;
//             }
//             _ASSERT(pos2 <= s_len);
//             m_CurrentGapLength += pos2 - pos;
//             m_CurrentGapChar = toupper(c);
//             pos = pos2 - 1; // `- 1` compensates for the `++pos` in the `for`
//             break;
//         }
//         case eCharType_JustIgnore:
//             break;
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2038-2043
// ```c++
//     NStr::TruncateSpacesInPlace(processed_title);
//     if (!processed_title.empty()) {
//         auto pDesc = Ref(new CSeqdesc());
//         pDesc->SetTitle() = processed_title;
//         bioseq.SetDescr().Set().push_back(std::move(pDesc));
//     }
// ```
pub fn parse_fasta(data: &[u8], protein: bool, lowercase: bool) -> Result<Vec<FastaRecord>> {
    parse_fasta_batches(data, protein, lowercase, usize::MAX, &mut |_| Ok(()))
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140-171
// ```c++
//     while (size_read < GetBatchSize()) {
//
//         if (End())
//             break;
//
//         CRef<CBlastSearchQuery> q;
//         try { q.Reset(m_Source->GetNextSequence(scope)); }
//         catch (const CObjReaderParseException& e) {
//             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
//                 break;
//             }
//             throw;
//         }
//         catch (const exception&) {
//             continue; //SB-2307. ignore well formed, not found accession
//         }
//
//         CConstRef<CSeq_loc> loc = q->GetQuerySeqLoc();
//
//         if (loc->IsInt()) {
//             size_read += sequence::GetLength(loc->GetInt().GetId(),
//                                              q->GetScope());
//         } else if (loc->IsWhole()) {
//             size_read += sequence::GetLength(loc->GetWhole(), q->GetScope());
//         } else {
//             // programmer error, CBlastInputSource should only return Seq-locs
//             // of type interval or whole
//             abort();
//         }
//
//         retval->AddQuery(q);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:260-293
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
// ```
pub(crate) fn parse_fasta_batches(
    data: &[u8],
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    parse_fasta_batches_with_warnings(data, protein, lowercase, batch_size, deliver, &mut |_| {
        Ok(())
    })
}
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2213-2217
// ```c++
//     if ( ! pMessageListener && (_eSeverity) <= eDiag_Warning ) {
//         LOG_POST_X(1, Warning << pLineExpt->Message());
//     } else if ( ! pMessageListener || ! pMessageListener->PutError( *pLineExpt ) )
//     {
//         throw CObjReaderParseException(DIAG_COMPILE_INFO, 0, _eErrCode, _MessageStrmOps, _uLineNum, _eSeverity);
// ```
// Native caller posts each reader event immediately; retained warnings continue
// to serve the independent B input-state probe. Batch delivery is unchanged.
pub(crate) fn parse_fasta_batches_with_warnings(
    data: &[u8],
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    parse_fasta_lines(
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:302-314
        // ```c++
        // CBlastFastaInputSource::CBlastFastaInputSource(const string& user_input,
        //                                        const CBlastInputSourceConfig& iconfig)
        //     : m_Config(iconfig),
        //       m_ReadProteins(iconfig.IsProteinInput())
        // {
        //     if (user_input.empty()) {
        //         NCBI_THROW(CInputException, eEmptyUserInput,
        //                    "No sequence input was provided");
        //     }
        //     m_LineReader.Reset(new CMemoryLineReader(user_input.c_str(),
        //                                              user_input.size()));
        //     x_InitInputReader();
        // }
        // ```
        FastaMemoryLines { remaining: data },
        protein,
        lowercase,
        batch_size,
        deliver,
        warn,
    )
}
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:350-375
// ```c++
//         if (c == '>' ) {
//             CTempString next_line = *++GetLineReader();
//             string strmodified;
//             if( NStr::StartsWith(next_line, ">?_") ) {
//                 CTempString tmp = next_line.substr(3);
//                 strmodified = ">";
//                 strmodified.append(tmp.data(), tmp.length());
//                 next_line = strmodified;
//             }
//             if( NStr::StartsWith(next_line, ">?") ) {
//                 // This is actually a data line. an assembly gap, in particular, which
//                 // we handle farther below
//                 GetLineReader().UngetLine();
//             } else {
//                 if(need_defline) {
//                     ParseDefLine(next_line, pMessageListener);
//                     need_defline = false;
//                     continue;
//                 } else {
//                     GetLineReader().UngetLine();
//                     // start of the next sequence
//                     break;
//                 }
//             }
//         }
//
// ```
fn parse_fasta_lines<L: AsRef<[u8]>>(
    lines: impl Iterator<Item = L>,
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    let mut batch_start = 0;
    let mut batch_length = 0;

    let mut records = Vec::new();
    let mut current: Option<FastaRecord> = None;
    let mut pending_title_warning = None;
    let mut raw_residues = 0;
    for (line_index, raw) in lines.enumerate() {
        let raw = raw.as_ref();
        let line = raw.strip_suffix(b"\r").unwrap_or(raw);
        // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:350-375
        // ```c++
        //         if (c == '>' ) {
        //             CTempString next_line = *++GetLineReader();
        //             string strmodified;
        //             if( NStr::StartsWith(next_line, ">?_") ) {
        //                 CTempString tmp = next_line.substr(3);
        //                 strmodified = ">";
        //                 strmodified.append(tmp.data(), tmp.length());
        //                 next_line = strmodified;
        //             }
        //             if( NStr::StartsWith(next_line, ">?") ) {
        //                 // This is actually a data line. an assembly gap, in particular, which
        //                 // we handle farther below
        //                 GetLineReader().UngetLine();
        //             } else {
        //                 if(need_defline) {
        //                     ParseDefLine(next_line, pMessageListener);
        //                     need_defline = false;
        //                     continue;
        //                 } else {
        //                     GetLineReader().UngetLine();
        //                     // start of the next sequence
        //                     break;
        //                 }
        //             }
        //         }
        //
        // ```
        let escaped_header = line
            .strip_prefix(b">?_")
            .map(|suffix| [b">".as_slice(), suffix].concat());
        let header = escaped_header.as_deref().unwrap_or(line);
        let is_header = line.starts_with(b">") && !header.starts_with(b">?");
        let line = if is_header {
            header
        } else {
            trim_input_space(line)
        };
        // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:399-408
        // ```c++
        //
        //         if ( !TestFlag(fNoSeqData) ) {
        //             try {
        //                 string strmodified;
        //                 if( NStr::StartsWith(line, ">?_") ) {
        //                     CTempString tmp = line.substr(3);
        //                     strmodified = ">";
        //                     strmodified.append(tmp.data(), tmp.length());
        //                     line = strmodified;
        //                 }
        // ```
        let rewritten_data;
        let line = if !is_header && line.starts_with(b">?_") {
            rewritten_data = [b">".as_slice(), &line[3..]].concat();
            rewritten_data.as_slice()
        } else {
            line
        };
        if line.is_empty() || matches!(line[0], b';' | b'#' | b'!') {
            continue;
        }
        if line.starts_with(b">?") {
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1098-1122
            // ```c++
            //
            //     // just in case there's a gap before this one,
            //     // even though somewhat unusual
            //     CloseGap();
            //
            //     // sRemainingLine will hold the part of the line left to parse
            //     TStr sRemainingLine = line.substr(2);
            //     NStr::TruncateSpacesInPlace(sRemainingLine);
            //
            //     const TSeqPos uPos = GetCurrentPos(eRawPos);
            //
            //     // check if size is unknown
            //     SGap::EKnownSize eIsKnown = SGap::eKnownSize_Yes;
            //     if( NStr::StartsWith(sRemainingLine, "unk") ) {
            //         eIsKnown = SGap::eKnownSize_No;
            //         sRemainingLine = sRemainingLine.substr(3);
            //         NStr::TruncateSpacesInPlace(sRemainingLine, NStr::eTrunc_Begin);
            //     }
            //
            //     // extract the gap size
            //     TSeqPos uGapSize = 0;
            //     {
            //         // find how many digits in number
            //         TStr::size_type uNumDigits = 0;
            //         while( uNumDigits != sRemainingLine.size() &&
            // ```
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:77-82
            // ```c++
            //     /// Override this method to force the parent class to ignore gaps
            //     /// @param len length of the gap? @sa CFastaReader
            //     protected:
            //     virtual void x_CloseGap(TSeqPos /*len*/, bool /*atStartOfLine*/,
            //                             ILineErrorListener * /*pMessageListener*/)
            //         { }
            // ```
            // Explicit assembly gaps become N/X in BLAST's default no-parse-gap mode.
            // Rich gap modifiers remain explicitly outside the local FASTA scope.
            let gap = trim_input_space(&line[2..]);
            let gap = trim_input_start(gap.strip_prefix(b"unk").unwrap_or(gap));
            if !gap.is_empty()
                && gap.iter().all(u8::is_ascii_digit)
                && std::str::from_utf8(gap)
                    .ok()
                    .and_then(|s| s.parse::<u32>().ok())
                    .is_some_and(|n| n > 0)
            {
                let record = current
                    .get_or_insert_with(|| new_record(records.len(), protein, String::new()));
                // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1387-1394
                // ```c++
                //         // Encountered >? lines; substitute runs of Ns or Xs as appropriate.
                //         string    new_data;
                //         char      gap_char(inst.IsAa() ? 'X' : 'N');
                //         SIZE_TYPE pos = 0;
                //         new_data.reserve(GetCurrentPos(ePosWithGaps));
                //         ITERATE (TGaps, it, m_Gaps) {
                //             // since we're not parsing gaps, we have to throw out
                //             // any specified extra information that can't be
                // ```
                // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1411-1425
                // ```c++
                //             if ((*it)->m_uPos > pos) {
                //                 new_data.append(m_SeqData, pos, (*it)->m_uPos - pos);
                //                 pos = (*it)->m_uPos;
                //             }
                //             new_data.append((*it)->m_uLen, gap_char);
                //         }
                //         if (m_CurrentPos > pos) {
                //             new_data.append(m_SeqData, pos, m_CurrentPos - pos);
                //         }
                //         swap(m_SeqData, new_data);
                //         m_Gaps.clear();
                //         m_CurrentPos += m_TotalGapLength;
                //         m_TotalGapLength = 0;
                //         m_CurrentGapChar = '\0';
                //     }
                // ```
                let count = std::str::from_utf8(gap).unwrap().parse::<u32>().unwrap() as usize;
                if count > i32::MAX as usize - record.sequence.len() {
                    bail!("BLASTX sequence exceeds Int4 length");
                }
                // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1082-1090
                // ```c++
                //     m_MaskRangeStart = GetCurrentPos(ePosWithGapsAndSegs);
                // }
                //
                // void CFastaReader::x_CloseMask(void)
                // {
                //     _ASSERT(m_MaskRangeStart != kInvalidSeqPos);
                //     m_CurrentMask->SetPacked_int().AddInterval
                //         (GetBestID(), m_MaskRangeStart, GetCurrentPos(ePosWithGapsAndSegs) - 1,
                //          eNa_strand_plus);
                // ```
                // NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:505-514
                // ```c++
                //     TSeqPos pos = m_CurrentPos;
                //     switch (pos_type) {
                //     case ePosWithGapsAndSegs:
                //         return pos + m_SegmentBase + m_TotalGapLength;
                //     case ePosWithGaps:
                //         return pos + m_TotalGapLength;
                //     case eRawPos:
                //         return pos;
                //     default:
                //         return kInvalidSeqPos;
                // ```
                if let Some(last) = record
                    .lowercase_masks
                    .last_mut()
                    .filter(|r| r.1 + 1 == record.sequence.len() as i32)
                {
                    last.1 += count as i32;
                }
                record.sequence.resize(
                    record.sequence.len() + count,
                    if protein { b'X' } else { b'N' },
                );
                continue;
            }
            bail!("unsupported FASTA assembly-gap size/modifiers: outside the declared local FASTA scope");
        }
        if is_header {
            if let Some(record) = current.take() {
                finish_record(record, &mut records, pending_title_warning.take())?;
                emit_fasta_batch(
                    &records,
                    batch_size,
                    &mut batch_start,
                    &mut batch_length,
                    deliver,
                )?;
            }
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:353-359
            // ```c++
            //             if( NStr::StartsWith(next_line, ">?_") ) {
            //                 CTempString tmp = next_line.substr(3);
            //                 strmodified = ">";
            //                 strmodified.append(tmp.data(), tmp.length());
            //                 next_line = strmodified;
            //             }
            //             if( NStr::StartsWith(next_line, ">?") ) {
            // ```
            let title_bytes = trim_input_start(&line[1..]);
            let end = title_bytes
                .iter()
                .skip(1)
                .position(|&b| b < b' ')
                .map_or(title_bytes.len(), |i| i + 1);
            let title = String::from_utf8(trim_input_end(&title_bytes[..end]).to_vec())
                .context("invalid UTF-8 FASTA title")?;
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1623-1677
            // ```c++
            //     // check for nuc or aa sequences at the end of the title
            //     const static size_t kWarnNumNucCharsAtEnd = 20;
            //     const static size_t kWarnAminoAcidCharsAtEnd = 50;
            //
            //     const size_t length = sLineText.length();
            //     SIZE_TYPE pos_to_check = length-1;
            //
            //     if((length > kWarnNumNucCharsAtEnd) && !TestFlag(fAssumeProt)) {
            //         // find last non-nuc character, within the last kWarnNumNucCharsAtEnd characters
            //         const SIZE_TYPE last_pos_to_check_for_nuc = (sLineText.length() - kWarnNumNucCharsAtEnd);
            //         for( ; pos_to_check >= last_pos_to_check_for_nuc; --pos_to_check ) {
            //             if( ! s_ASCII_IsUnAmbigNuc(sLineText[pos_to_check]) ) {
            //                 // found a character which is not an unambiguous nucleotide
            //                 break;
            //             }
            //         }
            //         if( pos_to_check < last_pos_to_check_for_nuc ) {
            //             FASTA_WARNING(iLineNum,
            //                 "FASTA-Reader: Title ends with at least " << kWarnNumNucCharsAtEnd
            //                 << " valid nucleotide characters.  Was the sequence "
            //                 << "accidentally put in the title line?",
            //                 ILineError::eProblem_UnexpectedNucResidues,
            //                 "defline"
            //                 );
            //             return true; // found problem
            //         }
            //     }
            //
            //     if((length > kWarnAminoAcidCharsAtEnd) && !TestFlag(fAssumeNuc)) {
            //         // check for aa's at the end of the title
            //         // for efficiency, continue where the nuc search left off, since
            //         // we know that nucs can be amino acids, also
            //         const SIZE_TYPE last_pos_to_check_for_amino_acid =
            //             ( sLineText.length() - kWarnAminoAcidCharsAtEnd );
            //         for( ; pos_to_check >= last_pos_to_check_for_amino_acid; --pos_to_check ) {
            //             // can't just use "isalpha" in case it includes characters
            //             // with diacritics (an accent, tilde, umlaut, etc.)
            //             const char ch = sLineText[pos_to_check];
            //             if( ( ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') ) {
            //                 // potential amino acid, so keep going
            //             } else {
            //                 // non-amino-acid found
            //                 break;
            //             }
            //         }
            //
            //         if( pos_to_check < last_pos_to_check_for_amino_acid ) {
            //             FASTA_WARNING(iLineNum,
            //                 "FASTA-Reader: Title ends with at least " << kWarnAminoAcidCharsAtEnd
            //                 << " valid amino acid characters.  Was the sequence "
            //                 << "accidentally put in the title line?",
            //                 ILineError::eProblem_UnexpectedAminoAcids,
            //                 "defline");
            //             return true; // found problem
            //         }
            // ```
            let raw_title = &title_bytes[..end];
            let kind = if !protein
                && raw_title.len() > 20
                && raw_title[raw_title.len() - 20..]
                    .iter()
                    .all(|c| b"ACGTacgt".contains(c))
            {
                Some("title_nucleotide_residues")
            } else if protein
                && raw_title.len() > 50
                && raw_title[raw_title.len() - 50..]
                    .iter()
                    .all(u8::is_ascii_alphabetic)
            {
                Some("title_amino_acids")
            } else {
                None
            };
            pending_title_warning = kind.map(|kind| InputWarning {
                line: line_index + 1,
                kind,
                positions: Vec::new(),
            });
            if let Some(warning) = &pending_title_warning {
                warn(warning)?;
            }
            current = Some(new_record(records.len(), protein, title));
            raw_residues = 0;
            continue;
        }
        let record =
            current.get_or_insert_with(|| new_record(records.len(), protein, String::new()));
        if raw_residues == 0 {
            check_data_line(line, line_index + 1)?;
        }
        let mut bad = Vec::new();
        let mut hyphens = false;
        for (pos, &c) in line.iter().enumerate() {
            if c == b';' {
                break;
            }
            if input_space(c) {
                continue;
            }
            if c == b'-' {
                hyphens = true;
                continue;
            }
            let upper = c.to_ascii_uppercase();
            let allowed = if protein {
                b"ABCDEFGHIJKLMNOPQRSTUVWXYZ*".contains(&upper)
            } else {
                b"ABCDGHKMNRSTUVWY".contains(&upper)
            };
            if !allowed {
                bad.push(pos + 1);
                continue;
            }
            raw_residues += 1;
            let offset = record.sequence.len() as i32;
            record.sequence.push(if !protein && upper == b'U' {
                b'T'
            } else {
                upper
            });
            if lowercase && c.is_ascii_lowercase() {
                if let Some(last) = record
                    .lowercase_masks
                    .last_mut()
                    .filter(|r| r.1 + 1 == offset)
                {
                    last.1 = offset;
                } else {
                    record.lowercase_masks.push((offset, offset));
                }
            }
        }
        if hyphens {
            let warning = InputWarning {
                line: line_index + 1,
                kind: "hyphens",
                positions: Vec::new(),
            };
            warn(&warning)?;
            record.warnings.push(warning);
        }
        if !bad.is_empty() {
            let warning = InputWarning {
                line: line_index + 1,
                kind: "invalid_residues",
                positions: bad,
            };
            warn(&warning)?;
            record.warnings.push(warning);
        }
    }
    if let Some(record) = current {
        finish_record(record, &mut records, pending_title_warning.take())?;
        emit_fasta_batch(
            &records,
            batch_size,
            &mut batch_start,
            &mut batch_length,
            deliver,
        )?;
    }
    if batch_start < records.len() {
        deliver(&records[batch_start..])?;
    }
    Ok(records)
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:364-367
// ```c++
//     CRef<CSeqIdGenerator> idgen
//         (new CSeqIdGenerator(m_Config.GetLocalIdCounterInitValue(),
//                              m_Config.GetLocalIdPrefix()));
//     m_InputReader->SetIDGenerator(*idgen);
// ```
fn new_record(index: usize, protein: bool, title: String) -> FastaRecord {
    FastaRecord {
        internal_id: format!(
            "{}_{n}",
            if protein { "Subject" } else { "Query" },
            n = index + 1
        ),
        title,
        sequence: Vec::new(),
        lowercase_masks: Vec::new(),
        warnings: Vec::new(),
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:437-447
// ```c++
//         ? 0 : m_Config.GetRange().GetTo();
//
//     // Get the sequence length
//     const TSeqPos seqlen = seq_entry->GetSeq().GetInst().GetLength();
//     //if (seqlen == 0) {
//     //    NCBI_THROW(CInputException, eEmptyUserInput,
//     //               "Query contains no sequence data");
//     //}
//     _ASSERT(seqlen != numeric_limits<TSeqPos>::max());
//     if (to > 0 && to < from) {
//         NCBI_THROW(CInputException, eInvalidRange,
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:457-461
// ```c++
//
//     // set sequence range
//     retval->SetInt().SetFrom(from);
//     retval->SetInt().SetTo((to > 0 && to < seqlen) ? to : (seqlen-1));
//
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1361-1365
// ```c++
//
//     // apply source mods *after* figuring out mol type
//     ITERATE(vector<SLineTextAndLoc>, title_ci, m_CurrentSeqTitles) {
//         ParseTitle(*title_ci, pMessageListener);
//     }
// ```
fn finish_record(
    mut record: FastaRecord,
    records: &mut Vec<FastaRecord>,
    title_warning: Option<InputWarning>,
) -> Result<()> {
    if let Some(warning) = title_warning {
        record.warnings.push(warning);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:334-340
    // ```c++
    //     }
    //     // This is necessary to enable the ignoring of gaps in classes derived from
    //     // CFastaReader
    //
    //    	flags+= CFastaReader::fHyphensIgnoreAndWarn;
    //
    //     flags+= CFastaReader::fDisableNoResidues;
    // ```
    // Empty records remain input objects; query/subject setup owns their
    if record.sequence.len() > i32::MAX as usize {
        bail!("BLASTX sequence exceeds Int4 length");
    }
    records.push(record);
    Ok(())
}
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:726-763
// ```c++
//         ( TestFlag(fAssumeNuc) && TestFlag(fForceType) ) ||
//         ( m_CurrentSeq && m_CurrentSeq->IsSetInst() &&
//           m_CurrentSeq->GetInst().IsSetMol() &&  m_CurrentSeq->IsNa() ) );
//     size_t ambig_nuc = 0;
//     for (size_t pos = 0;  pos < len_to_check;  ++pos) {
//         unsigned char c = s[pos];
//         if (s_ASCII_IsAlpha(c) ||  c == '*') {
//             ++good;
//             if( bIsNuc && s_ASCII_IsAmbigNuc(c) ) {
//                 ++ambig_nuc;
//             }
//         } else if( c == '-' ) {
//             if( ! bIgnoreHyphens ) {
//                 ++good;
//             }
//             // if bIgnoreHyphens == true, the "hyphens are ignored" warning
//             // will be triggered elsewhere
//         } else if (isspace(c)  ||  (c >= '0' && c <= '9')) {
//             // treat whitespace and digits as neutral
//         } else if (c == ';') {
//             break; // comment -- ignore rest of line
//         } else {
//             ++bad;
//         }
//     }
//     if (bad >= good / 3  &&
//         (len_to_check > 3  ||  good == 0  ||  bad > good))
//     {
//         FASTA_ERROR( LineNumber(),
//             "CFastaReader: Near line " << LineNumber()
//             << ", there's a line that doesn't look like plausible data, "
//             "but it's not marked as defline or comment.",
//             CObjReaderParseException::eFormat);
//     }
//     // warn if more than a certain percentage is ambiguous nucleotides
//     const static size_t kWarnPercentAmbiguous = 40; // e.g. "40" means "40%"
//     const size_t percent_ambig = (good == 0)?100:((ambig_nuc * 100) / good);
//     if( len_to_check > 3 && percent_ambig > kWarnPercentAmbiguous ) {
// ```
fn check_data_line(line: &[u8], number: usize) -> Result<()> {
    let len = line.len().min(70);
    let (mut good, mut bad) = (0, 0);
    for &c in &line[..len] {
        if c.is_ascii_alphabetic() || c == b'*' {
            good += 1;
        } else if c == b';' {
            break;
        } else if c != b'-' && !input_space(c) && !c.is_ascii_digit() {
            bad += 1;
        }
    }
    if bad >= good / 3 && (len > 3 || good == 0 || bad > good) {
        return Err(FastaParseError(format!("CFastaReader: Near line {number}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.")).into());
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/include/objtools/readers/reader_exception.hpp:46-49
// ```c++
// class CObjReaderParseException :
//     public CParseTemplException<CObjReaderException>
// {
// public:
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:754-758
// ```c++
//         FASTA_ERROR( LineNumber(),
//             "CFastaReader: Near line " << LineNumber()
//             << ", there's a line that doesn't look like plausible data, "
//             "but it's not marked as defline or comment.",
//             CObjReaderParseException::eFormat);
// ```
#[derive(Debug)]
pub(crate) struct FastaParseError(pub String);
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:181-184
// ```c++
//     catch (const CObjReaderParseException& e) {                             \
//         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
impl std::fmt::Display for FastaParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
// NCBI reference (598d8ae6): c++/include/objtools/readers/reader_exception.hpp:46-49
// ```c++
// class CObjReaderParseException :
//     public CParseTemplException<CObjReaderException>
// {
// public:
// ```
impl std::error::Error for FastaParseError {}

#[cfg(test)]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:86-94
// ```c++
//         for (i = 0; i < retval->last_context + 1; i++) {
//             retval->contexts[i].query_index =
//                 Blast_GetQueryIndexFromContext(i, program);
//             ASSERT(retval->contexts[i].query_index != -1);
//
//             retval->contexts[i].frame = BLAST_ContextToFrame(program,  i);
//             ASSERT(retval->contexts[i].frame != INT1_MAX);
//
//             retval->contexts[i].is_valid = TRUE;
// ```
mod tests {
    use super::*;
    #[test]
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:736-751
    // ```c++
    //             }
    //         } else if( c == '-' ) {
    //             if( ! bIgnoreHyphens ) {
    //                 ++good;
    //             }
    //             // if bIgnoreHyphens == true, the "hyphens are ignored" warning
    //             // will be triggered elsewhere
    //         } else if (isspace(c)  ||  (c >= '0' && c <= '9')) {
    //             // treat whitespace and digits as neutral
    //         } else if (c == ';') {
    //             break; // comment -- ignore rest of line
    //         } else {
    //             ++bad;
    //         }
    //     }
    //     if (bad >= good / 3  &&
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:974-990
    // ```c++
    //         case eCharType_Bad:
    //             if( bad_pos_line_num < 0 ) {
    //                 bad_pos_line_num = LineNumber();
    //             }
    //             bad_pos_vec.push_back(pos);
    //             break;
    //         default:
    //             _TROUBLE;
    //         }
    //     }
    //
    //     m_SeqData.resize(m_CurrentPos);
    //
    //     if( bIgnorableHyphenSeen ) {
    //         _ASSERT( bHyphensIgnoreAndWarn );
    //         FASTA_WARNING_EX(LineNumber(),
    //             "CFastaReader: Hyphens are invalid and will be ignored around line " << LineNumber(),
    // ```
    fn ncbi_reader_preserves_normalized_positions_and_warnings() {
        let r = parse_fasta(b">  q description   \r\natg1-2-3ACT;rest\r\n", false, true).unwrap();
        assert_eq!(r[0].title, "q description");
        assert_eq!(r[0].internal_id, "Query_1");
        assert_eq!(r[0].sequence, b"ATGACT");
        assert_eq!(r[0].lowercase_masks, vec![(0, 2)]);
        assert_eq!(r[0].warnings.len(), 2);
        assert_eq!(r[0].warnings[1].positions, vec![4, 6, 8]);
        assert!(parse_fasta(b">q\nATG?!ACT", false, false).is_err());
    }
    #[test]
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:382-390
    // ```c++
    //
    //         if (c == '!'  ||  c == '#' || c == ';') {
    //             // no content, just a comment or blank line
    //             continue;
    //         } else if (need_defline) {
    //             if (TestFlag(fDLOptional)) {
    //                 ParseDefLine(">", pMessageListener);
    //                 need_defline = false;
    //             } else {
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:437-442
    // ```c++
    //         ? 0 : m_Config.GetRange().GetTo();
    //
    //     // Get the sequence length
    //     const TSeqPos seqlen = seq_entry->GetSeq().GetInst().GetLength();
    //     //if (seqlen == 0) {
    //     //    NCBI_THROW(CInputException, eEmptyUserInput,
    // ```
    fn empty_file_record_and_optional_defline_differ() {
        assert!(parse_fasta(b"", false, false).unwrap().is_empty());
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:339-340
        // ```c++
        //
        //     flags+= CFastaReader::fDisableNoResidues;
        // ```
        assert!(parse_fasta(b">q\n", false, false).unwrap()[0]
            .sequence
            .is_empty());
        assert_eq!(
            parse_fasta(b"ATGACT", false, false).unwrap()[0].sequence,
            b"ATGACT"
        );
        let r = parse_fasta(b">dup\nATG\n>dup\nACT", false, false).unwrap();
        assert_eq!(
            r.iter().map(|r| r.internal_id.as_str()).collect::<Vec<_>>(),
            vec!["Query_1", "Query_2"]
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140-171
// ```c++
//     while (size_read < GetBatchSize()) {
//
//         if (End())
//             break;
//
//         CRef<CBlastSearchQuery> q;
//         try { q.Reset(m_Source->GetNextSequence(scope)); }
//         catch (const CObjReaderParseException& e) {
//             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
//                 break;
//             }
//             throw;
//         }
//         catch (const exception&) {
//             continue; //SB-2307. ignore well formed, not found accession
//         }
//
//         CConstRef<CSeq_loc> loc = q->GetQuerySeqLoc();
//
//         if (loc->IsInt()) {
//             size_read += sequence::GetLength(loc->GetInt().GetId(),
//                                              q->GetScope());
//         } else if (loc->IsWhole()) {
//             size_read += sequence::GetLength(loc->GetWhole(), q->GetScope());
//         } else {
//             // programmer error, CBlastInputSource should only return Seq-locs
//             // of type interval or whole
//             abort();
//         }
//
//         retval->AddQuery(q);
//     }
// ```
fn emit_fasta_batch(
    records: &[FastaRecord],
    batch_size: usize,
    start: &mut usize,
    length: &mut usize,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
) -> Result<()> {
    *length += records
        .last()
        .expect("completed FASTA record")
        .sequence
        .len();
    if *length >= batch_size {
        deliver(&records[*start..])?;
        *start = records.len();
        *length = 0;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:310-313
// ```c++
//     }
//     m_LineReader.Reset(new CMemoryLineReader(user_input.c_str(),
//                                              user_input.size()));
//     x_InitInputReader();
// ```
struct FastaMemoryLines<'a> {
    remaining: &'a [u8],
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:353-381
// ```c++
// CMemoryLineReader& CMemoryLineReader::operator++(void)
// {
//     /* If at EOF - noop */
//     if (AtEOF()) {
//         m_Line = CTempString(NULL);
//         return *this;
//     }
//     const char* p = m_Pos;
//     if ( p == m_Line.begin() ) {
//         /* If after UngetLine(), line is already in buffer, so end is known*/
//         p = m_Line.end();
//     } else {
//         /* Line is in stream, go char by char until delimiters */
//         while ( p < m_End  &&  *p != '\r'  && *p != '\n' ) {
//             ++p;
//         }
//         m_Line = CTempString(m_Pos, p - m_Pos);
//     }
//     // skip over delimiters until the beginning of the next string
//     if (p + 1 < m_End  &&  *p == '\r'  &&  p[1] == '\n') {
//         m_Pos = p + 2;
//     } else if (p < m_End) {
//         m_Pos = p + 1;
//     } else { // no final line break
//         m_Pos = p;
//     }
//     ++m_LineNumber;
//     return *this;
// }
// ```
impl<'a> Iterator for FastaMemoryLines<'a> {
    type Item = &'a [u8];
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:353-381
    // ```c++
    // CMemoryLineReader& CMemoryLineReader::operator++(void)
    // {
    //     /* If at EOF - noop */
    //     if (AtEOF()) {
    //         m_Line = CTempString(NULL);
    //         return *this;
    //     }
    //     const char* p = m_Pos;
    //     if ( p == m_Line.begin() ) {
    //         /* If after UngetLine(), line is already in buffer, so end is known*/
    //         p = m_Line.end();
    //     } else {
    //         /* Line is in stream, go char by char until delimiters */
    //         while ( p < m_End  &&  *p != '\r'  && *p != '\n' ) {
    //             ++p;
    //         }
    //         m_Line = CTempString(m_Pos, p - m_Pos);
    //     }
    //     // skip over delimiters until the beginning of the next string
    //     if (p + 1 < m_End  &&  *p == '\r'  &&  p[1] == '\n') {
    //         m_Pos = p + 2;
    //     } else if (p < m_End) {
    //         m_Pos = p + 1;
    //     } else { // no final line break
    //         m_Pos = p;
    //     }
    //     ++m_LineNumber;
    //     return *this;
    // }
    // ```
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining.is_empty() {
            return None;
        }
        let end = self
            .remaining
            .iter()
            .position(|&c| c == b'\r' || c == b'\n')
            .unwrap_or(self.remaining.len());
        let line = &self.remaining[..end];
        let next = if self.remaining.get(end..end + 2) == Some(b"\r\n") {
            end + 2
        } else if end < self.remaining.len() {
            end + 1
        } else {
            end
        };
        self.remaining = &self.remaining[next..];
        Some(line)
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
pub(crate) fn stream_is_empty<R: Read + Seek>(stream: &mut FastaStream<R>) -> bool {
    let input = &mut stream.input;
    let position = match input.stream_position() {
        Ok(position) => position,
        Err(_) => return false,
    };
    let mut byte = [0];
    loop {
        match input.read(&mut byte) {
            Ok(0) | Err(_) => return true,
            Ok(_) if input_space(byte[0]) => {}
            Ok(_) => {
                // IsIStreamEmpty restores position only after successful extraction.
                // NCBI unconditionally clears seekg failure and restores orig_state.
                let _ = input.seek(SeekFrom::Start(position));
                return false;
            }
        }
    }
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:166-174
// ```c++
//
//     switch (m_EOLStyle) {
//     case eEOL_unknown: x_AdvanceEOLUnknown();                   break;
//     case eEOL_cr:      x_AdvanceEOLSimple('\r', '\n');          break;
//     case eEOL_lf:      x_AdvanceEOLSimple('\n', '\r');          break;
//     case eEOL_crlf:    x_AdvanceEOLCRLF();                      break;
//     case eEOL_mixed:   NcbiGetline(*m_Stream, m_Line, "\r\n");  break;
//     }
//     return *this;
// ```
#[derive(Clone, Copy, PartialEq)]
enum EolStyle {
    Unknown,
    Cr,
    Lf,
    CrLf,
    Mixed,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
// ```c++
//       m_UngetLine(false), m_AutoEOL(eol_style == eEOL_unknown),
//       m_EOLStyle(eol_style)
// {
// }
//
//
// CStreamLineReader::CStreamLineReader(CNcbiIstream& is,
//                                      EOwnership ownership)
//     : m_Stream(&is, ownership), m_LineNumber(0), m_LastReadSize(0),
// ```
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:437-443
// ```c++
//     if (!del_ptr  &&  how != ePushback_NoCopy) {
//         del_ptr = new CT_CHAR_TYPE[buf_size];
//         buf = (CT_CHAR_TYPE*) memcpy(del_ptr, buf, buf_size);
//     }
//
//     (void) new CPushback_Streambuf(is, buf, buf_size, del_ptr);
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:134-143
// ```c++
// CPushback_Streambuf::CPushback_Streambuf(istream&      is,
//                                          CT_CHAR_TYPE* buf,
//                                          size_t        buf_size,
//                                          void*         del_ptr)
//     : m_Is(is), m_Next(0), m_Buf(buf), m_BufSize(buf_size), m_DelPtr(del_ptr)
// {
//     _ASSERT(m_Buf  &&  m_BufSize);
//     setp(0, 0);  // unbuffered output at this level of streambuf's hierarchy
//     setg(m_Buf, m_Buf, m_Buf + m_BufSize);
//     m_Sb = m_Is.rdbuf(this);
// ```
struct PushbackBuffer {
    bytes: Vec<u8>,
    position: usize,
    capacity: usize,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
// ```c++
//       m_UngetLine(false), m_AutoEOL(eol_style == eEOL_unknown),
//       m_EOLStyle(eol_style)
// {
// }
//
//
// CStreamLineReader::CStreamLineReader(CNcbiIstream& is,
//                                      EOwnership ownership)
//     : m_Stream(&is, ownership), m_LineNumber(0), m_LastReadSize(0),
// ```
pub(crate) struct FastaStream<R: Read> {
    input: BufReader<R>,
    failed: bool,
    eol: EolStyle,
    pushback: Vec<PushbackBuffer>,
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
    // ```c++
    //     x_FillBuffer((size_t) m_Sb->in_avail());
    //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
    // ```
    backend_available: Option<fn(&mut R) -> std::io::Result<usize>>,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
// ```c++
//       m_UngetLine(false), m_AutoEOL(eol_style == eEOL_unknown),
//       m_EOLStyle(eol_style)
// {
// }
//
//
// CStreamLineReader::CStreamLineReader(CNcbiIstream& is,
//                                      EOwnership ownership)
//     : m_Stream(&is, ownership), m_LineNumber(0), m_LastReadSize(0),
// ```
impl<R: Read> FastaStream<R> {
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
    // ```c++
    //       m_UngetLine(false), m_AutoEOL(eol_style == eEOL_unknown),
    //       m_EOLStyle(eol_style)
    // {
    // }
    //
    //
    // CStreamLineReader::CStreamLineReader(CNcbiIstream& is,
    //                                      EOwnership ownership)
    //     : m_Stream(&is, ownership), m_LineNumber(0), m_LastReadSize(0),
    // ```
    pub(crate) fn new(input: R) -> Self {
        Self {
            // NCBI reference (598d8ae6): c++/include/corelib/ncbistre.hpp:437-440
            // ```c++
            // #else
            // /// Portable alias for ifstream.
            // typedef IO_PREFIX::ifstream      CNcbiIfstream;
            // #endif
            // ```
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
            // ```c++
            // CT_INT_TYPE CPushback_Streambuf::underflow(void)
            // {
            //     // we are here because there is no more data in the pushback buffer
            //     _ASSERT(gptr()  &&  gptr() >= egptr());
            //
            // #ifdef NCBI_COMPILER_MIPSPRO
            //     if (m_MIPSPRO_ReadsomeGptrSetLevel  &&  m_MIPSPRO_ReadsomeGptr != gptr())
            //         return CT_EOF;
            //     m_MIPSPRO_ReadsomeGptr = (CT_CHAR_TYPE*)(-1L);
            // #endif //NCBI_COMPILER_MIPSPRO
            //
            //     x_FillBuffer((size_t) m_Sb->in_avail());
            //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
            // ```
            // The pinned CNcbiIfstream exposes 8191 readable bytes per default
            // window (independent C++ file_window oracle; tested at EOF/pushback
            // boundaries). Preserve that window because reuse retains EOF state.
            input: BufReader::with_capacity(8191, input),
            failed: false,
            eol: EolStyle::Unknown,
            pushback: Vec::new(),
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
            // ```c++
            //     x_FillBuffer((size_t) m_Sb->in_avail());
            //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
            // ```
            backend_available: None,
        }
    }
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
    // ```c++
    // CT_INT_TYPE CPushback_Streambuf::underflow(void)
    // {
    //     // we are here because there is no more data in the pushback buffer
    //     _ASSERT(gptr()  &&  gptr() >= egptr());
    //
    // #ifdef NCBI_COMPILER_MIPSPRO
    //     if (m_MIPSPRO_ReadsomeGptrSetLevel  &&  m_MIPSPRO_ReadsomeGptr != gptr())
    //         return CT_EOF;
    //     m_MIPSPRO_ReadsomeGptr = (CT_CHAR_TYPE*)(-1L);
    // #endif //NCBI_COMPILER_MIPSPRO
    //
    //     x_FillBuffer((size_t) m_Sb->in_avail());
    //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:306-358
    // ```c++
    // void CPushback_Streambuf::x_FillBuffer(size_t max_size)
    // {
    //     _ASSERT(m_Sb);
    //     if ( !max_size ) {
    //         ++max_size;
    //     }
    //
    //     CPushback_Streambuf* sb = dynamic_cast<CPushback_Streambuf*> (m_Sb);
    //     if ( !sb ) {
    //         CT_CHAR_TYPE* bp = 0;
    //         size_t buf_size = m_DelPtr
    //             ? (size_t)(m_Buf - (CT_CHAR_TYPE*) m_DelPtr) + m_BufSize : 0;
    //         if (buf_size < kMinBufSize) {
    //             buf_size = kMinBufSize;
    //             bp = new CT_CHAR_TYPE[buf_size];
    //         }
    //         streamsize r = (streamsize)(buf_size < max_size ? buf_size : max_size);
    //         streamsize n = m_Sb->sgetn(bp ? bp : (CT_CHAR_TYPE*) m_DelPtr, r);
    //         if (n <= 0) {
    //             // NB: For unknown reasons WorkShop6 can return -1 from sgetn :-/
    //             delete[] bp;
    //             return;
    //         }
    //         if (bp) {
    //             delete[] (CT_CHAR_TYPE*) m_DelPtr;
    //             m_DelPtr = bp;
    //         }
    //         m_Buf = (CT_CHAR_TYPE*) m_DelPtr;
    //         m_BufSize = buf_size;
    //         setg(m_Buf, m_Buf, m_Buf + n);
    //         return;
    //     }
    //
    //     _ASSERT(&m_Is  == &sb->m_Is);
    //     _ASSERT(m_Next == sb);
    //     m_Sb       = sb->m_Sb;
    //     m_Next     = sb->m_Next;
    //     sb->m_Sb   = 0;
    //     sb->m_Next = 0;
    //     if (sb->gptr() >= sb->egptr()) {
    //         delete sb;
    //         x_FillBuffer(max_size);
    //         return;
    //     }
    //     delete[] (CT_CHAR_TYPE*) m_DelPtr;
    //     m_Buf        = sb->m_Buf;
    //     m_BufSize    = sb->m_BufSize;
    //     m_DelPtr     = sb->m_DelPtr;
    //     sb->m_DelPtr = 0;
    //     setg(sb->gptr(), sb->gptr(), sb->egptr());
    //     delete sb;
    // }
    //
    // ```
    fn raw_peek(&mut self) -> Option<u8> {
        while let Some(buffer) = self.pushback.last() {
            if let Some(&byte) = buffer.bytes.get(buffer.position) {
                return Some(byte);
            }
            if self.pushback.len() > 1 {
                let below = self.pushback.remove(self.pushback.len() - 2);
                if below.position < below.bytes.len() {
                    *self.pushback.last_mut().unwrap() = below;
                }
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:306-338
            // ```c++
            // void CPushback_Streambuf::x_FillBuffer(size_t max_size)
            // {
            //     _ASSERT(m_Sb);
            //     if ( !max_size ) {
            //         ++max_size;
            //     }
            //
            //     CPushback_Streambuf* sb = dynamic_cast<CPushback_Streambuf*> (m_Sb);
            //     if ( !sb ) {
            //         CT_CHAR_TYPE* bp = 0;
            //         size_t buf_size = m_DelPtr
            //             ? (size_t)(m_Buf - (CT_CHAR_TYPE*) m_DelPtr) + m_BufSize : 0;
            //         if (buf_size < kMinBufSize) {
            //             buf_size = kMinBufSize;
            //             bp = new CT_CHAR_TYPE[buf_size];
            //         }
            //         streamsize r = (streamsize)(buf_size < max_size ? buf_size : max_size);
            //         streamsize n = m_Sb->sgetn(bp ? bp : (CT_CHAR_TYPE*) m_DelPtr, r);
            //         if (n <= 0) {
            //             // NB: For unknown reasons WorkShop6 can return -1 from sgetn :-/
            //             delete[] bp;
            //             return;
            //         }
            //         if (bp) {
            //             delete[] (CT_CHAR_TYPE*) m_DelPtr;
            //             m_DelPtr = bp;
            //         }
            //         m_Buf = (CT_CHAR_TYPE*) m_DelPtr;
            //         m_BufSize = buf_size;
            //         setg(m_Buf, m_Buf, m_Buf + n);
            //         return;
            //     }
            //
            // ```
            let buffered = self.input.buffer().len();
            let available = if buffered != 0 {
                buffered
            } else if let Some(in_avail) = self.backend_available {
                in_avail(self.input.get_mut()).unwrap_or(0)
            } else {
                0
            };
            let capacity = self.pushback.last().unwrap().capacity.max(4096);
            let requested = capacity.min(available.max(1));
            let mut bytes = vec![0; requested];
            let mut length = 0;
            while length < requested {
                match self.input.read(&mut bytes[length..]) {
                    Ok(0) | Err(_) => break,
                    Ok(n) => length += n,
                }
            }
            if length == 0 {
                return None;
            }
            bytes.truncate(length);
            let buffer = self.pushback.last_mut().unwrap();
            buffer.capacity = capacity;
            buffer.bytes = bytes;
            buffer.position = 0;
        }
        match self.input.fill_buf() {
            Ok(bytes) => bytes.first().copied(),
            Err(_) => None,
        }
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:87-98
    // ```c++
    //             iostate = NcbiEofbit;
    //             break;
    //         }
    //         SIZE_TYPE delim_pos = delims.find(CT_TO_CHAR_TYPE(ch));
    //         if (delim_pos != NPOS) {
    //             // Special case -- if two different delimiters are back to
    //             // back and in the same order as in delims, treat them as
    //             // a single delimiter (necessary for correct handling of
    //             // DOS/MAC-style CR/LF endings).
    //             ch = is.rdbuf()->sgetc();
    //             if (!CT_EQ_INT_TYPE(ch, CT_EOF)
    //                 &&  delims.find(CT_TO_CHAR_TYPE(ch), delim_pos + 1) != NPOS) {
    // ```
    fn raw_take(&mut self) -> Option<u8> {
        let byte = self.raw_peek()?;
        if let Some(buffer) = self.pushback.last_mut() {
            buffer.position += 1;
        } else {
            self.input.consume(1);
        }
        Some(byte)
    }
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:413-443
    // ```c++
    //         if (how == ePushback_Stepback
    //             ||  (how == ePushback_Copy
    //                  &&  buf_size <= (del_ptr
    //                                   ? CPushback_Streambuf::kMinBufSize
    //                                   : CPushback_Streambuf::kMinBufSize >> 4))) {
    //             CT_CHAR_TYPE* bp = sb->gptr();
    //             size_t avail = bp - sb->m_Buf;
    //             size_t take  = avail < buf_size ? avail : buf_size;
    //             if (take) {
    //                 bp -= take;
    //                 buf_size -= take;
    //                 if (how != ePushback_Stepback  &&  bp != buf + buf_size) {
    //                     memmove(bp, buf + buf_size, take);
    //                 }
    //                 sb->setg(bp, bp, sb->egptr());
    //             }
    //         }
    //     }
    //
    //     if ( !buf_size ) {
    //         delete[] (CT_CHAR_TYPE*) del_ptr;
    //         return;
    //     }
    //
    //     if (!del_ptr  &&  how != ePushback_NoCopy) {
    //         del_ptr = new CT_CHAR_TYPE[buf_size];
    //         buf = (CT_CHAR_TYPE*) memcpy(del_ptr, buf, buf_size);
    //     }
    //
    //     (void) new CPushback_Streambuf(is, buf, buf_size, del_ptr);
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:134-143
    // ```c++
    // CPushback_Streambuf::CPushback_Streambuf(istream&      is,
    //                                          CT_CHAR_TYPE* buf,
    //                                          size_t        buf_size,
    //                                          void*         del_ptr)
    //     : m_Is(is), m_Next(0), m_Buf(buf), m_BufSize(buf_size), m_DelPtr(del_ptr)
    // {
    //     _ASSERT(m_Buf  &&  m_BufSize);
    //     setp(0, 0);  // unbuffered output at this level of streambuf's hierarchy
    //     setg(m_Buf, m_Buf, m_Buf + m_BufSize);
    //     m_Sb = m_Is.rdbuf(this);
    // ```
    fn push_back(&mut self, bytes: &[u8]) {
        let mut remaining = bytes.len();
        if let Some(buffer) = self.pushback.last_mut() {
            if bytes.len() <= (4096 >> 4) {
                let take = buffer.position.min(remaining);
                buffer.position -= take;
                remaining -= take;
                buffer.bytes[buffer.position..buffer.position + take]
                    .copy_from_slice(&bytes[remaining..]);
            }
        }
        if remaining != 0 {
            self.pushback.push(PushbackBuffer {
                bytes: bytes[..remaining].to_vec(),
                position: 0,
                capacity: remaining,
            });
            // Installing a new streambuf clears state; recycling the old one does not.
            self.failed = false;
        }
    }
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    // ```c++
    // bool CStreamLineReader::AtEOF(void) const
    // {
    //     return !m_UngetLine &&
    //         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
    // }
    // ```
    pub(crate) fn at_eof(&mut self) -> bool {
        if self.failed {
            return true;
        }
        if self.raw_peek().is_none() {
            self.failed = true;
        }
        self.failed
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:151-166
    // ```c++
    //     SIZE_TYPE size = 0;
    //     SIZE_TYPE max_size = str.max_size();
    //     do {
    //         CT_INT_TYPE nextc = is.get();
    //         if (CT_EQ_INT_TYPE(nextc, CT_EOF)
    //             ||  CT_EQ_INT_TYPE(nextc, CT_TO_INT_TYPE(delim))) {
    //             ++size;
    //             break;
    //         }
    //         if ( !is.unget() )
    //             break;
    //         if (size == max_size) {
    //             is.clear(NcbiFailbit);
    //             break;
    //         }
    //         SIZE_TYPE n = max_size - size;
    // ```
    fn next_byte(&mut self) -> Option<u8> {
        if self.failed {
            return None;
        }
        let byte = self.raw_take();
        if byte.is_none() {
            self.failed = true;
        }
        byte
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:87-105
    // ```c++
    //             iostate = NcbiEofbit;
    //             break;
    //         }
    //         SIZE_TYPE delim_pos = delims.find(CT_TO_CHAR_TYPE(ch));
    //         if (delim_pos != NPOS) {
    //             // Special case -- if two different delimiters are back to
    //             // back and in the same order as in delims, treat them as
    //             // a single delimiter (necessary for correct handling of
    //             // DOS/MAC-style CR/LF endings).
    //             ch = is.rdbuf()->sgetc();
    //             if (!CT_EQ_INT_TYPE(ch, CT_EOF)
    //                 &&  delims.find(CT_TO_CHAR_TYPE(ch), delim_pos + 1) != NPOS) {
    //                 is.rdbuf()->sbumpc();
    //                 delim_count = 2;
    //             } else {
    //                 delim_count = 1;
    //             }
    //             break;
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:151-174
    // ```c++
    //     SIZE_TYPE size = 0;
    //     SIZE_TYPE max_size = str.max_size();
    //     do {
    //         CT_INT_TYPE nextc = is.get();
    //         if (CT_EQ_INT_TYPE(nextc, CT_EOF)
    //             ||  CT_EQ_INT_TYPE(nextc, CT_TO_INT_TYPE(delim))) {
    //             ++size;
    //             break;
    //         }
    //         if ( !is.unget() )
    //             break;
    //         if (size == max_size) {
    //             is.clear(NcbiFailbit);
    //             break;
    //         }
    //         SIZE_TYPE n = max_size - size;
    //         is.get(buf, n < sizeof(buf) ? n : sizeof(buf), delim);
    //         n = (size_t) is.gcount();
    //         str.append(buf, n);
    //         size += n;
    //         _ASSERT(size == str.length());
    //     } while ( is.good() );
    // #endif
    //
    // ```
    fn getline(&mut self, delimiters: &[u8]) -> (Vec<u8>, Option<u8>) {
        let mut line = Vec::new();
        let mut last = None;
        loop {
            let byte = if delimiters.len() == 1 {
                self.next_byte()
            } else {
                self.raw_take()
            };
            let Some(byte) = byte else {
                self.failed = true;
                break;
            };
            last = Some(byte);
            if let Some(position) = delimiters.iter().position(|&delimiter| delimiter == byte) {
                if delimiters.len() > 1 {
                    if let Some(next) = self.raw_peek() {
                        if delimiters[position + 1..].contains(&next) {
                            last = self.raw_take();
                        }
                    }
                }
                break;
            }
            line.push(byte);
        }
        (line, last)
    }
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:245-267
    // ```c++
    // CStreamLineReader::EEOLStyle CStreamLineReader::x_AdvanceEOLSimple(char eol,
    //                                                                    char alt_eol)
    // {
    //     SIZE_TYPE pos;
    //     NcbiGetline(*m_Stream, m_Line, eol, &m_LastReadSize);
    //     if (m_AutoEOL  &&  (pos = m_Line.find(alt_eol)) != NPOS) {
    //         ++pos;
    //         if (eol != '\n'  ||  pos != m_Line.size()) {
    //             // an *immediately* preceding CR is quite all right
    //             CStreamUtils::Pushback(*m_Stream, m_Line.data() + pos,
    //                                    m_Line.size() - pos);
    //             m_EOLStyle = eEOL_mixed;
    //         }
    //         m_Line.resize(pos - 1);
    //         m_LastReadSize = pos;
    //         return (m_EOLStyle == eEOL_mixed) ? m_EOLStyle : eEOL_crlf;
    //     } else if (m_AutoEOL  &&  eol == '\r'  &&
    //                CT_EQ_INT_TYPE(m_Stream->peek(), CT_TO_INT_TYPE(alt_eol))) {
    //         m_Stream->get();
    //         ++m_LastReadSize;
    //         return eEOL_crlf;
    //     }
    //     return (eol == '\r') ? eEOL_cr : eEOL_lf;
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:437-443
    // ```c++
    //     if (!del_ptr  &&  how != ePushback_NoCopy) {
    //         del_ptr = new CT_CHAR_TYPE[buf_size];
    //         buf = (CT_CHAR_TYPE*) memcpy(del_ptr, buf, buf_size);
    //     }
    //
    //     (void) new CPushback_Streambuf(is, buf, buf_size, del_ptr);
    // }
    // ```
    fn advance_simple(&mut self, eol: u8, alternate: u8) -> (Vec<u8>, EolStyle) {
        let (mut line, _) = self.getline(&[eol]);
        if let Some(position) = line.iter().position(|&byte| byte == alternate) {
            let position = position + 1;
            if eol != b'\n' || position != line.len() {
                self.push_back(&line[position..]);
                self.eol = EolStyle::Mixed;
            }
            line.truncate(position - 1);
            let style = if self.eol == EolStyle::Mixed {
                EolStyle::Mixed
            } else {
                EolStyle::CrLf
            };
            return (line, style);
        }
        if eol == b'\r' && !self.failed && self.raw_peek() == Some(alternate) {
            self.raw_take();
            return (line, EolStyle::CrLf);
        }
        (
            line,
            if eol == b'\r' {
                EolStyle::Cr
            } else {
                EolStyle::Lf
            },
        )
    }
}
// NCBI reference (598d8ae6): c++/include/corelib/ncbistre.hpp:437-440
// ```c++
// #else
// /// Portable alias for ifstream.
// typedef IO_PREFIX::ifstream      CNcbiIfstream;
// #endif
// ```
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
// ```c++
// CT_INT_TYPE CPushback_Streambuf::underflow(void)
// {
//     // we are here because there is no more data in the pushback buffer
//     _ASSERT(gptr()  &&  gptr() >= egptr());
//
// #ifdef NCBI_COMPILER_MIPSPRO
//     if (m_MIPSPRO_ReadsomeGptrSetLevel  &&  m_MIPSPRO_ReadsomeGptr != gptr())
//         return CT_EOF;
//     m_MIPSPRO_ReadsomeGptr = (CT_CHAR_TYPE*)(-1L);
// #endif //NCBI_COMPILER_MIPSPRO
//
//     x_FillBuffer((size_t) m_Sb->in_avail());
//     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
// ```
// The registered std::ifstream backend reports buffered availability first,
// then remaining regular-file bytes without forcing the next 8191-byte window.
// Independent instrumented pinned CPushback_Streambuf calibrates this hint;
// sgetn may span windows when the preserved allocation is larger than 8191.
fn file_available(file: &mut std::fs::File) -> std::io::Result<usize> {
    let metadata = file.metadata()?;
    if metadata.is_file() {
        let position = file.stream_position()?;
        Ok(metadata.len().saturating_sub(position) as usize)
    } else {
        // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
        // ```c++
        //     x_FillBuffer((size_t) m_Sb->in_avail());
        //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
        // ```
        // The registered native CNcbiIfstream exposes queued pipe bytes via
        // FIONREAD before forcing a file-buffer refill (independent C++ oracle).
        #[cfg(unix)]
        {
            rustix::io::ioctl_fionread(file)
                .map(|n| n as usize)
                .map_err(Into::into)
        }
        #[cfg(not(unix))]
        {
            Ok(0)
        }
    }
}
// NCBI reference (598d8ae6): c++/include/corelib/ncbistre.hpp:437-440
// ```c++
// #else
// /// Portable alias for ifstream.
// typedef IO_PREFIX::ifstream      CNcbiIfstream;
// #endif
// ```
impl FastaStream<std::fs::File> {
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
    // ```c++
    // CT_INT_TYPE CPushback_Streambuf::underflow(void)
    // {
    //     // we are here because there is no more data in the pushback buffer
    //     _ASSERT(gptr()  &&  gptr() >= egptr());
    //
    // #ifdef NCBI_COMPILER_MIPSPRO
    //     if (m_MIPSPRO_ReadsomeGptrSetLevel  &&  m_MIPSPRO_ReadsomeGptr != gptr())
    //         return CT_EOF;
    //     m_MIPSPRO_ReadsomeGptr = (CT_CHAR_TYPE*)(-1L);
    // #endif //NCBI_COMPILER_MIPSPRO
    //
    //     x_FillBuffer((size_t) m_Sb->in_avail());
    //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
    // ```
    pub(crate) fn from_file(input: std::fs::File) -> Self {
        let mut stream = Self::new(input);
        stream.backend_available = Some(file_available);
        stream
    }
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:154-174
// ```c++
// CStreamLineReader& CStreamLineReader::operator++(void)
// {
//     /* If at EOF - noop */
//     if (AtEOF()) {
//         m_Line = string();
//         return *this;
//     }
//     ++m_LineNumber;
//     if ( m_UngetLine ) {
//         m_UngetLine = false;
//         return *this;
//     }
//
//     switch (m_EOLStyle) {
//     case eEOL_unknown: x_AdvanceEOLUnknown();                   break;
//     case eEOL_cr:      x_AdvanceEOLSimple('\r', '\n');          break;
//     case eEOL_lf:      x_AdvanceEOLSimple('\n', '\r');          break;
//     case eEOL_crlf:    x_AdvanceEOLCRLF();                      break;
//     case eEOL_mixed:   NcbiGetline(*m_Stream, m_Line, "\r\n");  break;
//     }
//     return *this;
// ```
impl<R: Read> Iterator for FastaStream<R> {
    type Item = Vec<u8>;
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:219-241
    // ```c++
    // CStreamLineReader::EEOLStyle CStreamLineReader::x_AdvanceEOLUnknown(void)
    // {
    //     _ASSERT(m_AutoEOL);
    //     NcbiGetline(*m_Stream, m_Line, "\r\n", &m_LastReadSize);
    //     m_Stream->unget();
    //     CT_INT_TYPE eol = m_Stream->get();
    //     if (CT_EQ_INT_TYPE(eol, CT_TO_INT_TYPE('\r'))) {
    //         m_EOLStyle = eEOL_cr;
    //     } else if (CT_EQ_INT_TYPE(eol, CT_TO_INT_TYPE('\n'))) {
    //         // NcbiGetline doesn't yield enough information to determine
    //         // whether eEOL_lf or eEOL_crlf is more appropriate, and not
    //         // all streams allow tellg() (which could otherwise resolve
    //         // matters), so defer further analysis to x_AdvanceEOLCRLF,
    //         // which will be responsible for reading the next line and
    //         // supports switching to eEOL_lf as appropriate.
    //         //
    //         // An alternative approach would have been to pass \n\r rather
    //         // than \r\n, and then check for an immediately following \n
    //         // if eol turned out to be \r, but that would miscount an
    //         // actual(!) \n\r sequence as a single line break.
    //         m_EOLStyle = eEOL_crlf;
    //     }
    //     return m_EOLStyle;
    // ```
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:271-281
    // ```c++
    // CStreamLineReader::EEOLStyle CStreamLineReader::x_AdvanceEOLCRLF(void)
    // {
    //     if (m_AutoEOL) {
    //         EEOLStyle style = x_AdvanceEOLSimple('\n', '\r');
    //         if (style == eEOL_mixed) {
    //             // found an embedded CR
    //             m_EOLStyle = eEOL_cr;
    //         } else if (style != eEOL_crlf) {
    //             m_EOLStyle = eEOL_lf;
    //         }
    //     } else {
    // ```
    fn next(&mut self) -> Option<Self::Item> {
        if self.at_eof() {
            return None;
        }
        let line = match self.eol {
            EolStyle::Unknown => {
                let (line, last) = self.getline(b"\r\n");
                // Successful unget clears eofbit before get rereads the last byte.
                if self.failed && !line.is_empty() {
                    self.failed = false;
                }
                match last {
                    Some(b'\r') => self.eol = EolStyle::Cr,
                    Some(b'\n') => self.eol = EolStyle::CrLf,
                    _ => {}
                }
                line
            }
            EolStyle::Cr => self.advance_simple(b'\r', b'\n').0,
            EolStyle::Lf => self.advance_simple(b'\n', b'\r').0,
            EolStyle::CrLf => {
                let (line, style) = self.advance_simple(b'\n', b'\r');
                if style == EolStyle::Mixed {
                    self.eol = EolStyle::Cr;
                } else if style != EolStyle::CrLf {
                    self.eol = EolStyle::Lf;
                }
                line
            }
            EolStyle::Mixed => self.getline(b"\r\n").0,
        };
        Some(line)
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-262
// ```c++
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
// 	    BLAST_PROF_START( APP.LOOP.PRE );
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140-150
// ```c++
//     while (size_read < GetBatchSize()) {
//
//         if (End())
//             break;
//
//         CRef<CBlastSearchQuery> q;
//         try { q.Reset(m_Source->GetNextSequence(scope)); }
//         catch (const CObjReaderParseException& e) {
//             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
//                 break;
//             }
// ```
pub(crate) fn parse_fasta_stream_with_warnings<R: Read>(
    stream: &mut FastaStream<R>,
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    if stream.at_eof() {
        return Ok(Vec::new());
    }
    let records = parse_fasta_lines(stream, protein, lowercase, batch_size, deliver, warn)?;
    if records.is_empty() {
        deliver(&[])?;
    }
    Ok(records)
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:930-931
// ```c++
//         case '\t': case '\n': case '\v': case '\f': case '\r': case ' ':
//             continue;
// ```
fn input_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3153-3161
// ```c++
//     SIZE_TYPE beg = 0;
//     if (where == NStr::eTrunc_Begin  ||  where == NStr::eTrunc_Both) {
//         _ASSERT(beg < length);
//         while ( isspace((unsigned char) str[beg]) ) {
//             if (++beg == length) {
//                 return empty_str;
//             }
//         }
//     }
// ```
fn trim_input_start(bytes: &[u8]) -> &[u8] {
    let start = bytes
        .iter()
        .position(|&byte| !input_space(byte))
        .unwrap_or(bytes.len());
    &bytes[start..]
}
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3162-3172
// ```c++
//     SIZE_TYPE end = length;
//     if ( where == NStr::eTrunc_End  ||  where == NStr::eTrunc_Both ) {
//         _ASSERT(beg < end);
//         while (isspace((unsigned char) str[--end])) {
//             if (beg == end) {
//                 return empty_str;
//             }
//         }
//         _ASSERT(beg <= end  &&  !isspace((unsigned char) str[end]));
//         ++end;
//     }
// ```
fn trim_input_end(bytes: &[u8]) -> &[u8] {
    let end = bytes
        .iter()
        .rposition(|&byte| !input_space(byte))
        .map_or(0, |p| p + 1);
    &bytes[..end]
}
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3187-3190
// ```c++
// CTempString NStr::TruncateSpaces_Unsafe(const CTempString str, ETrunc where)
// {
//     return s_TruncateSpaces(str, where, CTempString());
// }
// ```
fn trim_input_space(bytes: &[u8]) -> &[u8] {
    trim_input_end(trim_input_start(bytes))
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
// ```c++
// bool CStreamLineReader::AtEOF(void) const
// {
//     return !m_UngetLine &&
//         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
// }
// ```
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
#[cfg(test)]
mod stream_tests {
    use super::*;
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
    // ```c++
    //     x_FillBuffer((size_t) m_Sb->in_avail());
    //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
    // ```
    use std::io::Cursor;
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    // ```c++
    // bool CStreamLineReader::AtEOF(void) const
    // {
    //     return !m_UngetLine &&
    //         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
    // }
    // ```
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
    trait NativeRead: Read + Seek {
        // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
        // ```c++
        //     x_FillBuffer((size_t) m_Sb->in_avail());
        //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
        // ```
        fn available(&mut self) -> std::io::Result<usize>;
    }
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
    impl NativeRead for std::fs::File {
        // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
        // ```c++
        //     x_FillBuffer((size_t) m_Sb->in_avail());
        //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
        // ```
        fn available(&mut self) -> std::io::Result<usize> {
            file_available(self)
        }
    }
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    // ```c++
    // bool CStreamLineReader::AtEOF(void) const
    // {
    //     return !m_UngetLine &&
    //         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
    // }
    // ```
    struct Input {
        data: Cursor<Vec<u8>>,
        seekable: bool,
        fault: Option<u64>,
        error: bool,
    }
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    // ```c++
    // bool CStreamLineReader::AtEOF(void) const
    // {
    //     return !m_UngetLine &&
    //         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
    // }
    // ```
    impl Read for Input {
        // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
        // ```c++
        // bool CStreamLineReader::AtEOF(void) const
        // {
        //     return !m_UngetLine &&
        //         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
        // }
        // ```
        fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
            if self.fault.is_some_and(|p| self.data.position() >= p) {
                return if self.error {
                    Err(std::io::Error::other("filebuf read failure"))
                } else {
                    Ok(0)
                };
            }
            let length = bytes.len().min(1);
            self.data.read(&mut bytes[..length])
        }
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-871
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
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
    // ```c++
    //     x_FillBuffer((size_t) m_Sb->in_avail());
    //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
    // ```
    impl NativeRead for Input {
        // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
        // ```c++
        //     x_FillBuffer((size_t) m_Sb->in_avail());
        //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
        // ```
        fn available(&mut self) -> std::io::Result<usize> {
            Ok(0)
        }
    }
    impl Seek for Input {
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-871
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
        // ```
        fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> {
            if self.seekable {
                self.data.seek(position)
            } else {
                Err(std::io::Error::new(
                    std::io::ErrorKind::Unsupported,
                    "nonseekable stream",
                ))
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:154-174
    // ```c++
    // CStreamLineReader& CStreamLineReader::operator++(void)
    // {
    //     /* If at EOF - noop */
    //     if (AtEOF()) {
    //         m_Line = string();
    //         return *this;
    //     }
    //     ++m_LineNumber;
    //     if ( m_UngetLine ) {
    //         m_UngetLine = false;
    //         return *this;
    //     }
    //
    //     switch (m_EOLStyle) {
    //     case eEOL_unknown: x_AdvanceEOLUnknown();                   break;
    //     case eEOL_cr:      x_AdvanceEOLSimple('\r', '\n');          break;
    //     case eEOL_lf:      x_AdvanceEOLSimple('\n', '\r');          break;
    //     case eEOL_crlf:    x_AdvanceEOLCRLF();                      break;
    //     case eEOL_mixed:   NcbiGetline(*m_Stream, m_Line, "\r\n");  break;
    //     }
    //     return *this;
    // ```
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
    #[test]
    fn pinned_cpp_stream_state_and_line_oracle() {
        // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:353-381
        // ```c++
        // CMemoryLineReader& CMemoryLineReader::operator++(void)
        // {
        //     /* If at EOF - noop */
        //     if (AtEOF()) {
        //         m_Line = CTempString(NULL);
        //         return *this;
        //     }
        //     const char* p = m_Pos;
        //     if ( p == m_Line.begin() ) {
        //         /* If after UngetLine(), line is already in buffer, so end is known*/
        //         p = m_Line.end();
        //     } else {
        //         /* Line is in stream, go char by char until delimiters */
        //         while ( p < m_End  &&  *p != '\r'  && *p != '\n' ) {
        //             ++p;
        //         }
        //         m_Line = CTempString(m_Pos, p - m_Pos);
        //     }
        //     // skip over delimiters until the beginning of the next string
        //     if (p + 1 < m_End  &&  *p == '\r'  &&  p[1] == '\n') {
        //         m_Pos = p + 2;
        //     } else if (p < m_End) {
        //         m_Pos = p + 1;
        //     } else { // no final line break
        //         m_Pos = p;
        //     }
        //     ++m_LineNumber;
        //     return *this;
        // }
        // ```
        // Decode lossless independent C++ fixture bytes; no Rust-generated
        // expected lines or normalization of actual public report streams.
        let decode = |encoded: &str| -> Vec<u8> {
            if let Some(runs) = encoded.strip_prefix('~') {
                let mut bytes = Vec::new();
                for run in runs.split(',') {
                    let (byte, count) = run.split_once(':').unwrap();
                    let byte = u8::from_str_radix(byte, 16).unwrap();
                    bytes.extend(std::iter::repeat_n(byte, count.parse::<usize>().unwrap()));
                }
                bytes
            } else {
                encoded
                    .as_bytes()
                    .chunks(2)
                    .map(|b| u8::from_str_radix(std::str::from_utf8(b).unwrap(), 16).unwrap())
                    .collect()
            }
        };
        for row in include_str!("../../../tests/unit/blastx_stage_e_io_stream_expected.tsv").lines()
        {
            let fields: Vec<_> = row.split('\t').collect();
            assert_eq!(fields.len(), 9);
            let data = decode(fields[1]);
            let count = fields[8].parse::<usize>().unwrap();
            let expected: Vec<Vec<u8>> = if count == 0 {
                Vec::new()
            } else {
                fields[7].split('|').map(decode).collect()
            };
            assert_eq!(expected.len(), count);
            // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:353-381
            // ```c++
            // CMemoryLineReader& CMemoryLineReader::operator++(void)
            // {
            //     /* If at EOF - noop */
            //     if (AtEOF()) {
            //         m_Line = CTempString(NULL);
            //         return *this;
            //     }
            //     const char* p = m_Pos;
            //     if ( p == m_Line.begin() ) {
            //         /* If after UngetLine(), line is already in buffer, so end is known*/
            //         p = m_Line.end();
            //     } else {
            //         /* Line is in stream, go char by char until delimiters */
            //         while ( p < m_End  &&  *p != '\r'  && *p != '\n' ) {
            //             ++p;
            //         }
            //         m_Line = CTempString(m_Pos, p - m_Pos);
            //     }
            //     // skip over delimiters until the beginning of the next string
            //     if (p + 1 < m_End  &&  *p == '\r'  &&  p[1] == '\n') {
            //         m_Pos = p + 2;
            //     } else if (p < m_End) {
            //         m_Pos = p + 1;
            //     } else { // no final line break
            //         m_Pos = p;
            //     }
            //     ++m_LineNumber;
            //     return *this;
            // }
            // ```
            if fields[0].starts_with("memory_") {
                let lines: Vec<_> = FastaMemoryLines { remaining: &data }.collect();
                assert_eq!(
                    lines,
                    expected.iter().map(Vec::as_slice).collect::<Vec<_>>(),
                    "{} memory lines",
                    fields[0]
                );
                continue;
            }
            let fault: i64 = fields[4].parse().unwrap();
            let file_backend = fields[0].starts_with("file_");
            for &error in if file_backend {
                &[false][..]
            } else {
                &[false, true][..]
            } {
                let temporary = std::env::temp_dir().join(format!(
                    "losatx-e-stream-unit-{}-{}",
                    std::process::id(),
                    fields[0]
                ));
                let input: Box<dyn NativeRead> = if file_backend {
                    assert!(!temporary.exists());
                    std::fs::write(&temporary, &data).unwrap();
                    Box::new(std::fs::File::open(&temporary).unwrap())
                } else {
                    Box::new(Input {
                        data: Cursor::new(data.clone()),
                        seekable: fields[2] == "1",
                        fault: if fault < 0 { None } else { Some(fault as u64) },
                        error,
                    })
                };
                let mut stream = FastaStream::new(input);
                // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
                // ```c++
                //     x_FillBuffer((size_t) m_Sb->in_avail());
                //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
                // ```
                stream.backend_available = Some(|input| input.available());
                let empty = fields[3] == "1" && stream_is_empty(&mut stream);
                assert_eq!(empty, fields[6] == "1", "{} empty error={error}", fields[0]);
                let lines: Vec<Vec<u8>> = if empty { Vec::new() } else { stream.collect() };
                if file_backend {
                    std::fs::remove_file(&temporary).unwrap();
                }
                assert_eq!(lines, expected, "{} lines error={error}", fields[0]);
            }
        }
    }
}
