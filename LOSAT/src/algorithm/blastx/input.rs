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

use crate::blastinput::fasta_reader::RecordSource;
use crate::blastinput::seq_range::{
    record_interval, RecordInterval, RecordPlacement, SequenceRange,
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
    /// The defline line as read, for the BLAST archive (outfmt 11); `None` (8 bytes, no
    /// allocation) unless the reader was asked to keep it. A record without a defline line
    /// has none either (the archive then writes `>` and the title).
    pub source: Option<Box<crate::blastinput::fasta_reader::RecordSource>>,
}

impl FastaRecord {
    /// The Seq-ids `-parse_deflines` read, `None` without the option
    /// (`fasta_reader::FastaRecord::believed`).
    pub fn believed(&self) -> Option<&crate::blastinput::fasta_reader::believed::BelievedIds> {
        self.source.as_ref()?.believed.as_deref()
    }

    /// What `Query=` and `# Query:` show (`fasta_reader::FastaRecord::query_line`).
    pub fn query_line(&self) -> std::borrow::Cow<'_, str> {
        if self.believed().is_none() || self.title.is_empty() {
            return match self.believed() {
                Some(_) => std::borrow::Cow::Borrowed(&self.internal_id),
                None => std::borrow::Cow::Borrowed(&self.title),
            };
        }
        std::borrow::Cow::Owned(format!("{} {}", self.internal_id, self.title))
    }

    /// The id of the record in NCBI's warnings (`fasta_reader::FastaRecord::warning_id`).
    pub fn warning_id(&self) -> &[u8] {
        match self.believed() {
            Some(ids) => &ids.warning_id,
            None => self.internal_id.as_bytes(),
        }
    }

    /// The Seq-id of the record in the Seq-aligns and the archive.
    pub fn seq_id(&self) -> crate::report::seqid::SeqId<'_> {
        match self.believed() {
            Some(ids) => crate::report::seqid::SeqId::Db(&ids.best),
            None => crate::report::seqid::SeqId::Local(&self.internal_id),
        }
    }

    /// Where the letters lie in the record: the start of the `-query_loc` / `-subject_loc`
    /// interval that `cut_in_place` cut the record to, and the record's length (0 and the
    /// letters' length for a record searched whole).
    pub fn placement(&self) -> RecordPlacement {
        match self.source.as_deref() {
            Some(RecordSource {
                cut_from,
                uncut_length,
                ..
            }) if *uncut_length > 0 => RecordPlacement {
                offset: *cut_from,
                length: *uncut_length,
            },
            _ => RecordPlacement {
                offset: 0,
                length: self.sequence.len(),
            },
        }
    }

    /// The letters of the record as the archive's Bioseq holds them: the whole record when a range
    /// cut it keeping them (`cut_in_place`'s `keep_whole`), else the searched letters.
    pub fn archive_letters(&self) -> &[u8] {
        self.source
            .as_deref()
            .and_then(|source| source.uncut.as_deref())
            .unwrap_or(&self.sequence)
    }

    /// The record's lower-case runs `[left, right]` as the reader found them: in the whole record
    /// for a record that keeps its whole letters (those of `archive_letters`), else in the searched
    /// letters (`lowercase_masks`).
    pub fn archive_lowercase_runs(&self) -> Vec<(i32, i32)> {
        match self
            .source
            .as_deref()
            .and_then(|source| source.uncut.as_deref())
        {
            Some(uncut) => crate::algorithm::common::masks::lowercase_runs(uncut),
            None => self.lowercase_masks.clone(),
        }
    }

    /// Whether the record has lower-case runs (`archive_lowercase_runs`).
    pub fn has_lowercase_runs(&self) -> bool {
        match self
            .source
            .as_deref()
            .and_then(|source| source.uncut.as_deref())
        {
            Some(uncut) => uncut.iter().any(u8::is_ascii_lowercase),
            None => !self.lowercase_masks.is_empty(),
        }
    }

    /// The record as the Seq-loc of the archive's `seq-loc-list` (`-parse_deflines`, an id that is
    /// not local): the searched interval (`fasta_reader::FastaRecord::archive_seq_loc`).
    pub fn archive_seq_loc(&self) -> Option<crate::report::archive::bioseq::QuerySeqLoc<'_>> {
        let from = self.placement().offset as i32;
        self.believed()?
            .seq_loc(from, from + self.sequence.len() as i32 - 1)
    }

    /// Whether the letter before `from` lies in a lower-case run. NCBI restricts the runs of the
    /// whole record to each chunk of a split query (`split::x_chunk_masks`): of the letters
    /// outside an interval that starts at `from`, only that one can reach a chunk's masks (the
    /// first chunk's first letter), so a record cut there keeps its whole letters only then.
    pub fn lower_case_before(&self, from: usize) -> bool {
        let Some(letter) = from.checked_sub(1) else {
            return false;
        };
        self.lowercase_masks
            .iter()
            .any(|&(left, right)| left as usize <= letter && letter <= right as usize)
    }

    /// Cuts the record to its letters `[from, to_exclusive)` in place, as NCBI searches a record
    /// with a range (`seq_range::record_cut`): the lower-case masks are clipped to the letters
    /// and moved with them, and the source keeps where they lie in the record (`placement`).
    /// `keep_whole`: the source keeps the record's whole letters too (`archive_letters`,
    /// `archive_lowercase_runs`), for the archive's Bioseq (outfmt 11), the lower case of the
    /// whole record that outfmt 0 shows on a `-subject_loc` subject's rows, and the run of the
    /// letter before the interval that the chunks of a split query restrict
    /// (`lower_case_before`, `split::x_chunk_masks`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1041-1046
    /// ```c
    ///    for (seqloc = *mask; seqloc; seqloc = next_loc) {
    ///       next_loc = seqloc->next;
    ///       seqloc->ssr->left = MAX(0, seqloc->ssr->left - from);
    ///       seqloc->ssr->right = MIN(seqloc->ssr->right, to) - from;
    ///       /* If this mask location does not intersect the [from,to] interval,
    ///          do not add it to the newly constructed list and free its contents. */
    /// ```
    pub fn cut_in_place(&mut self, from: usize, to_exclusive: usize, keep_whole: bool) {
        let length = self.sequence.len();
        let source = self.source.get_or_insert_with(Default::default);
        source.uncut_length = length;
        source.cut_from = from;
        if keep_whole {
            // The searched letters are copied out and the record's letters become the whole
            // letters, so the record holds them once. As the shared reader keeps them: lower case
            // where the lower-case masks lie, which the archive's `LCaseMask` reads back
            // (`archive_lowercase_runs`).
            let mut uncut = std::mem::take(&mut self.sequence);
            self.sequence = uncut[from..to_exclusive].to_vec();
            for &(left, right) in &self.lowercase_masks {
                uncut[left as usize..=right as usize].make_ascii_lowercase();
            }
            source.uncut = Some(uncut);
        } else {
            self.sequence.truncate(to_exclusive);
            self.sequence.drain(..from);
        }
        let (from, last) = (from as i32, to_exclusive as i32 - 1);
        self.lowercase_masks.retain_mut(|(left, right)| {
            *left = (*left - from).max(0);
            *right = (*right).min(last) - from;
            *left <= *right
        });
    }
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
    /// The message of kind `defline_id` (the Seq-id warnings of `-parse_deflines`) and of kind
    /// `gap_line` (the warnings of a `>?` gap line, `fasta_reader::reader::gap_line_size`).
    pub text: Vec<u8>,
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
        false,
        None,
        None,
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
    keep_defline: bool,
    parse_deflines: Option<BelievedRead<'_>>,
    range: Option<&SequenceRange>,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    let mut believed_keys = std::collections::HashSet::new();
    let mut believed_generated = 0;
    let mut believed_molecule = None;
    let mut batch = ReaderBatch {
        range,
        keep_whole: keep_defline,
        ..ReaderBatch::default()
    };

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
        // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:383-398
        // ```c++
        //         if (c == '!'  ||  c == '#' || c == ';') {
        //             // no content, just a comment or blank line
        //             continue;
        //         } else if (need_defline) {
        //             if (TestFlag(fDLOptional)) {
        //                 ParseDefLine(">", pMessageListener);
        //                 need_defline = false;
        //             } else {
        //                 const auto lineNum = LineNumber();
        //                 GetLineReader().UngetLine();
        //                 NCBI_THROW2(CObjReaderParseException, eNoDefline,
        //                             "CFastaReader: Input doesn't start with"
        //                             " a defline or comment around line " + NStr::NumericToString(lineNum),
        //                              lineNum);
        //             }
        //         }
        // ```
        // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:774-780
        // ```c++
        // void CFastaReader::ParseDataLine(
        //     const TStr& s, ILineErrorListener * pMessageListener)
        // {
        //     if( NStr::StartsWith(s, ">?") ) {
        //         ParseGapLine(s, pMessageListener);
        //         return;
        //     }
        // ```
        // `-parse_deflines` reads without `fDLOptional`, so a first line that is not a defline,
        // a comment or blank stops the read before `ParseDataLine`: a `>?N` gap line too, and
        // whatever its size or modifiers.
        if !is_header && parse_deflines.is_some() && current.is_none() {
            return Err(NoDeflineError(line_index as u64 + 1).into());
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
            // The gap line is read as the reader of the other programs reads it
            // (`fasta_reader::reader::gap_line_size`, the whole `ParseGapLine`): a size that is not
            // a number from 1 to 4294967295 is 1 after `Bad gap size`, and the modifiers only warn.
            // The warnings go out as the line is read, before those of the data lines after it.
            let mut messages = Vec::new();
            let count = crate::blastinput::fasta_reader::reader::gap_line_size(
                line,
                line_index as u64 + 1,
                &mut |message| {
                    messages.push(message.strip_suffix(b"\n").unwrap_or(message).to_vec());
                    Ok(())
                },
            )? as usize;
            let record = current.get_or_insert_with(|| {
                new_record(records.len() + batch.skipped, protein, String::new())
            });
            for text in messages {
                let warning = InputWarning {
                    line: line_index + 1,
                    kind: "gap_line",
                    positions: Vec::new(),
                    text,
                };
                warn(&warning)?;
                record.warnings.push(warning);
            }
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
        if is_header {
            if let Some(record) = current.take() {
                finish_record(
                    record,
                    &mut records,
                    pending_title_warning.take(),
                    protein,
                    parse_deflines.and_then(|read| read.id_scope),
                    warn,
                )?;
                emit_fasta_batch(&mut records, batch_size, &mut batch, deliver)?;
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
            let believed = match parse_deflines {
                Some(read) => Some(read_believed_defline(
                    line,
                    protein,
                    read,
                    keep_defline,
                    line_index as u64 + 1,
                    &mut believed_keys,
                    &mut believed_generated,
                    warn,
                )?),
                None => None,
            };
            let title_bytes = match &believed {
                Some(defline) => trim_input_start(&line[defline.title_search..]),
                None => trim_input_start(&line[1..]),
            };
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
                text: Vec::new(),
            });
            let mut record = new_record(records.len() + batch.skipped, protein, title);
            believed_molecule = None;
            if keep_defline || believed.is_some() {
                record.source = Some(Box::new(crate::blastinput::fasta_reader::RecordSource {
                    defline: if keep_defline {
                        header.to_vec()
                    } else {
                        Vec::new()
                    },
                    uncut: None,
                    believed: believed.as_ref().map(|defline| defline.ids.clone()),
                    cut_from: 0,
                    uncut_length: 0,
                }));
            }
            if let Some(defline) = believed {
                record.internal_id = defline.ids.bare_id();
                believed_molecule = defline.molecule;
            }
            current = Some(record);
            raw_residues = 0;
            continue;
        }
        let record = current.get_or_insert_with(|| {
            new_record(records.len() + batch.skipped, protein, String::new())
        });
        if raw_residues == 0 {
            check_data_line(line, line_index + 1)?;
        }
        if believed_molecule.is_some_and(|named| named != protein)
            && !crate::blastinput::fasta_reader::believed::reads_the_same_in_both(line)
        {
            anyhow::bail!(
                "the Seq-id of the record on line {} names the other molecule under \
                 -parse_deflines (NCBI BLAST+ reads its letters as that molecule); this is not \
                 supported by LOSAT's BLASTX",
                line_index + 1
            );
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
                // The letters `ABCDEFGHIJKLMNOPQRSTUVWXYZ*`, tested without a table in memory:
                // `contains` on a static table first steps byte by byte to an aligned address,
                // so its cost per letter depended on where the linker placed the table.
                matches!(upper, b'A'..=b'Z' | b'*')
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
                text: Vec::new(),
            };
            warn(&warning)?;
            record.warnings.push(warning);
        }
        if !bad.is_empty() {
            let warning = InputWarning {
                line: line_index + 1,
                kind: "invalid_residues",
                positions: bad,
                text: crate::blastinput::fasta_reader::believed::nuc_or_prot(believed_molecule)
                    .as_bytes()
                    .to_vec(),
            };
            warn(&warning)?;
            record.warnings.push(warning);
        }
    }
    if let Some(record) = current {
        finish_record(
            record,
            &mut records,
            pending_title_warning.take(),
            protein,
            parse_deflines.and_then(|read| read.id_scope),
            warn,
        )?;
        emit_fasta_batch(&mut records, batch_size, &mut batch, deliver)?;
    }
    // A last batch of records that the range skipped is delivered empty (NCBI's query factory
    // then refuses it).
    if batch.read > 0 {
        deliver(&records[batch.start..])?;
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
/// How `-parse_deflines` reads: the strand of a query's Seq-loc, and the search's ids of the other
/// role's records (`FastaStream::id_scope`).
#[derive(Clone, Copy)]
struct BelievedRead<'a> {
    strand: crate::report::archive::bioseq::SeqLocStrand,
    id_scope: Option<&'a crate::seqid::id_scope::IdScope>,
}

/// A defline under `-parse_deflines` (`fasta_reader::believed`).
struct BelievedLine {
    title_search: usize,
    ids: std::sync::Arc<crate::blastinput::fasta_reader::believed::BelievedIds>,
    molecule: Option<bool>,
}

/// `fasta_reader::reader::parse_believed_def_line` for BLASTX's reader.
fn read_believed_defline(
    line: &[u8],
    protein: bool,
    read: BelievedRead<'_>,
    keep_defline: bool,
    line_number: u64,
    keys: &mut std::collections::HashSet<u64>,
    generated: &mut usize,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<BelievedLine> {
    use crate::blastinput::fasta_reader::believed::{
        duplicate_id, generated_ids, named_molecule, parse_believed_defline,
        refuse_range_in_archive, BelievedIds,
    };
    let mut written = Ok(());
    let parsed = parse_believed_defline(line, protein, "blastx", line_number, &mut |message| {
        if written.is_ok() {
            written = warn(&InputWarning {
                line: line_number as usize,
                kind: "defline_id",
                positions: Vec::new(),
                text: message.to_vec(),
            });
        }
    });
    written?;
    let parsed = parsed.map_err(FastaParseError)?;
    let role = if protein { "subject" } else { "query" };
    if parsed.range && keep_defline {
        return Err(refuse_range_in_archive(role, "BLASTX"));
    }
    let molecule = named_molecule(&parsed.ids);
    let ids = if parsed.ids.is_empty() || parsed.range {
        *generated += 1;
        generated_ids(&format!(
            "{}_{generated}",
            if protein { "Subject" } else { "Query" }
        ))
    } else {
        parsed.ids
    };
    if let Some(refusal) = duplicate_id(&ids, keys, role, "BLASTX") {
        return Err(refusal);
    }
    Ok(BelievedLine {
        title_search: parsed.title_search,
        ids: std::sync::Arc::new(BelievedIds::new(ids, read.strand)),
        molecule,
    })
}

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
        source: None,
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
/// NCBI posts a record's title warning when it assembles the record, after the warnings of its
/// data lines; a record whose data stops the reading posts none. Under `-parse_deflines` the
/// record's ids then join the search's ids (`fasta_reader`'s `assemble_seq`, whose
/// `add_to_id_scope` follows `parse_title`), so LOSAT's refusal of a record that a record of the
/// other role names comes after the record's messages.
fn finish_record(
    mut record: FastaRecord,
    records: &mut Vec<FastaRecord>,
    title_warning: Option<InputWarning>,
    protein: bool,
    id_scope: Option<&crate::seqid::id_scope::IdScope>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<()> {
    if let Some(warning) = title_warning {
        warn(&warning)?;
        record.warnings.push(warning);
    }
    // A nucleotide query and a protein subject are never one record (`RecordMark`).
    if let (Some(scope), Some(believed)) = (id_scope, record.believed()) {
        let role = if protein { "subject" } else { "query" };
        scope.add_record(&believed.ids, None, role, "BLASTX")?;
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

/// The reader's `eNoDefline` at a line (1-based): under `-parse_deflines` the first line that
/// is neither blank, a comment nor a defline. Its text is the reader's (`CFastaReader: ...`, a
/// query's `BLAST query error`, `native::input_error`); the subjects end with `GetAllSeqs`'s
/// error instead (`native::subject_input_error`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:391-396
/// ```c++
///                 const auto lineNum = LineNumber();
///                 GetLineReader().UngetLine();
///                 NCBI_THROW2(CObjReaderParseException, eNoDefline,
///                             "CFastaReader: Input doesn't start with"
///                             " a defline or comment around line " + NStr::NumericToString(lineNum),
///                              lineNum);
/// ```
#[derive(Debug)]
pub(crate) struct NoDeflineError(pub u64);

impl std::fmt::Display for NoDeflineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CFastaReader: Input doesn't start with a defline or comment around line {}",
            self.0
        )
    }
}

impl std::error::Error for NoDeflineError {}

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
///
/// With a range (`-query_loc`) the batch counts the length of the whole record and cuts the record
/// to its interval; a record whose interval starts more than one letter past its end adds nothing
/// and is not kept (it keeps its number: `ReaderBatch::skipped`).
fn emit_fasta_batch(
    records: &mut Vec<FastaRecord>,
    batch_size: usize,
    batch: &mut ReaderBatch<'_>,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
) -> Result<()> {
    let record = records.last_mut().expect("completed FASTA record");
    let length = record.sequence.len();
    batch.read += 1;
    match batch.range.map(|range| record_interval(range, length)) {
        None => batch.length += length,
        Some(RecordInterval::Letters { from, to }) => {
            // The letter before the interval is the only one outside it that a split query's
            // chunks read (`split::x_chunk_masks`); a record keeps its whole letters for it.
            let keep_whole = batch.keep_whole || record.lower_case_before(from);
            record.cut_in_place(from, to + 1, keep_whole);
            batch.length += length;
        }
        Some(RecordInterval::Empty) => {
            record.cut_in_place(length, length, batch.keep_whole);
            batch.length += length;
        }
        Some(RecordInterval::PastEnd) => {
            records.pop();
            batch.skipped += 1;
        }
    }
    if batch.length >= batch_size {
        deliver(&records[batch.start..])?;
        batch.start = records.len();
        batch.length = 0;
        batch.read = 0;
    }
    Ok(())
}

/// The reader batch that `emit_fasta_batch` fills.
#[derive(Default)]
struct ReaderBatch<'r> {
    /// The range that cuts each record (`-query_loc`), if any.
    range: Option<&'r SequenceRange>,
    /// Whether a cut record keeps its whole letters (the archive, `FastaRecord::cut_in_place`).
    keep_whole: bool,
    /// The first kept record of the batch.
    start: usize,
    /// The letters the batch has read (`size_read`).
    length: usize,
    /// The records the batch has read, kept or not.
    read: usize,
    /// The records the range skipped so far: NCBI's reader numbers them (`Query_<n>`).
    skipped: usize,
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
    if stream.no_position {
        return false;
    }
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
    /// Whether the input has no stream position, as a pipe has none (`tellg() < 0` of
    /// `IsIStreamEmpty`; `input_files::has_stream_position`): then it is never empty.
    no_position: bool,
    /// `-parse_deflines` (`fasta_reader::believed`).
    pub(crate) parse_deflines: bool,
    /// `-strand`, read with `-parse_deflines`.
    pub(crate) strand: crate::report::archive::bioseq::SeqLocStrand,
    /// The search whose records this stream reads under `-parse_deflines`.
    pub(crate) id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    /// With NCBI's data loaders on (`ncbi_environment::ApplicationSettings::data_loaders`), the
    /// role ("query", "subject") named in LOSAT's refusal of a first line that NCBI tries as a
    /// Seq-id (`fasta_reader::seq_id::reject_seq_id_line`); `None` reads every first line as
    /// FASTA, as NCBI does with the data loaders off.
    pub(crate) seq_id_role: Option<&'static str>,
    /// `CStreamLineReader::UngetLine`: the line that the next read gives again (the Seq-id step
    /// of `parse_fasta_stream_keeping` read it and handed it back).
    unget_line: Option<Vec<u8>>,
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
            no_position: false,
            parse_deflines: false,
            id_scope: None,
            strand: crate::report::archive::bioseq::SeqLocStrand::Both,
            seq_id_role: None,
            unget_line: None,
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
        if self.unget_line.is_some() {
            return false;
        }
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
    // Windows reports a pipe or a console as a file whose position can be read
    // (`input_files::has_stream_position`); elsewhere `stream_position` below fails on one.
    if metadata.is_file()
        && (cfg!(not(windows)) || crate::blastinput::input_files::has_stream_position(file))
    {
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
        let no_position = !crate::blastinput::input_files::has_stream_position(&input);
        let mut stream = Self::new(input);
        stream.no_position = no_position;
        stream.backend_available = Some(file_available);
        stream
    }
}

/// How the search reads its query file: `-parse_deflines`, the search's ids and `-strand`
/// (`FastaStream`'s fields). The reads of the file before the search
/// (`blast_xml::prescan_queries`, `archive::prescan_query_titles`) read it the same way, so that
/// they see the search's records and titles and stop at the reader error that stops the search
/// (under `-parse_deflines`, a first line that is not a defline).
pub(crate) struct QueryStreamSettings {
    pub(crate) parse_deflines: bool,
    pub(crate) id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    pub(crate) strand: crate::report::archive::bioseq::SeqLocStrand,
    /// NCBI's data loaders (`FastaStream::seq_id_role`).
    pub(crate) data_loaders: bool,
}

impl QueryStreamSettings {
    /// The search's stream of the query file.
    pub(crate) fn stream(&self, file: std::fs::File) -> FastaStream<std::fs::File> {
        let mut stream = FastaStream::from_file(file);
        stream.parse_deflines = self.parse_deflines;
        stream.id_scope = self.id_scope.clone();
        stream.strand = self.strand;
        stream.seq_id_role = self.data_loaders.then_some("query");
        stream
    }

    /// A stream of the query file for a read before the search: its ids go to a copy of the
    /// search's (`IdScope::copy`), so that it stops at the record the search refuses (a query
    /// with the id of a subject) and the search reads its records into its ids as before.
    pub(crate) fn stream_before_search(&self, file: std::fs::File) -> FastaStream<std::fs::File> {
        let mut stream = self.stream(file);
        stream.id_scope = self
            .id_scope
            .as_deref()
            .map(|scope| std::sync::Arc::new(scope.copy()));
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
        if let Some(line) = self.unget_line.take() {
            return Some(line);
        }
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
    parse_fasta_stream_keeping(
        stream, protein, lowercase, batch_size, false, None, deliver, warn,
    )
}
/// `parse_fasta_stream_with_warnings`, keeping each record's defline line (`keep_defline`: the
/// BLAST archive, outfmt 11, writes it) and cutting each record to `range` (`-query_loc`,
/// `emit_fasta_batch`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn parse_fasta_stream_keeping<R: Read>(
    stream: &mut FastaStream<R>,
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    keep_defline: bool,
    range: Option<&SequenceRange>,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    if stream.at_eof() {
        return Ok(Vec::new());
    }
    let id_scope = stream.id_scope.clone();
    let parse_deflines = stream.parse_deflines.then_some(BelievedRead {
        strand: stream.strand,
        id_scope: id_scope.as_deref(),
    });
    // A file without records gives one empty batch; one whose records a range skipped gave it
    // already.
    let mut delivered = false;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:130-162
    // ```c++
    //     virtual CRef<CSeq_entry> ReadOneSeq(ILineErrorListener * pMessageListener) {
    //
    //         const string line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
    //         if ( !line.empty() && isalnum(line.data()[0]&0xff) ) {
    //             try {
    //                 CRef<CSeq_id> id(new CSeq_id(line, (CSeq_id::fParse_AnyRaw |
    // 							CSeq_id::fParse_ValidLocal)));
    // ...
    //         // If all fails, fall back to parent's implementation
    //         GetLineReader().UngetLine();
    //         return CFastaReader::ReadOneSeq(pMessageListener);
    //     }
    // ```
    // BLASTX reads its query and its subjects with the reader of the other programs
    // (`CBlastFastaInputSource`), whose every record first tries its first line as a Seq-id when
    // the data loaders are on (`seq_id_role`). Only the input's first line can be one: a record
    // ends before the next `>` line. LOSAT refuses such a line as the other programs' reader does
    // (`fasta_reader::reader::read_one_seq`), here, before `CFastaReader`'s read (and its
    // `-parse_deflines` stop at a first line that is not a defline), which gets the line back
    // (`UngetLine`).
    if let Some(role) = stream.seq_id_role {
        if let Some(line) = stream.next() {
            if let Some(refusal) =
                crate::blastinput::fasta_reader::seq_id::reject_seq_id_line(&line, "BLASTX", role)
            {
                return Err(refusal);
            }
            stream.unget_line = Some(line);
        }
    }
    let records = parse_fasta_lines(
        stream,
        protein,
        lowercase,
        batch_size,
        keep_defline,
        parse_deflines,
        range,
        &mut |records| {
            delivered = true;
            deliver(records)
        },
        warn,
    )?;
    if !delivered {
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

#[cfg(test)]
mod range_tests {
    use super::*;
    use std::io::Cursor;

    /// The batches a query file of `data` gives with `range` and `batch_size`: for each, the ids,
    /// letters and placements of its records; and the records the reader keeps.
    fn batches(
        data: &[u8],
        range: &SequenceRange,
        batch_size: usize,
    ) -> (
        Vec<Vec<(String, Vec<u8>, RecordPlacement)>>,
        Vec<FastaRecord>,
    ) {
        let mut delivered = Vec::new();
        let kept = parse_fasta_stream_keeping(
            &mut FastaStream::new(Cursor::new(data.to_vec())),
            false,
            true,
            batch_size,
            false,
            Some(range),
            &mut |records| {
                delivered.push(
                    records
                        .iter()
                        .map(|record| {
                            (
                                record.internal_id.clone(),
                                record.sequence.clone(),
                                record.placement(),
                            )
                        })
                        .collect(),
                );
                Ok(())
            },
            &mut |_| Ok(()),
        )
        .unwrap();
        (delivered, kept)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:144-168
    // ```c++
    //         catch (const exception&) {
    //             continue; //SB-2307. ignore well formed, not found accession
    //         }
    //
    //         CConstRef<CSeq_loc> loc = q->GetQuerySeqLoc();
    //
    //         if (loc->IsInt()) {
    //             size_read += sequence::GetLength(loc->GetInt().GetId(),
    //                                              q->GetScope());
    // ```
    // A batch counts the whole length of each searched record; a skipped record (its interval
    // starts more than one letter past its end) counts nothing, is not delivered and keeps its
    // number; the letters delivered are the interval's.
    #[test]
    fn a_range_cuts_the_records_and_skips_the_short_ones() {
        let range = SequenceRange { from: 4, to: 7 };
        let data = b">a\nACGTACGTAC\n>b\nACG\n>c\nTTTTGGGG\n>d\nACGT\n>e\nACGTA\n";
        let (delivered, kept) = batches(data, &range, 15);
        let place = |offset, length| RecordPlacement { offset, length };
        assert_eq!(
            delivered,
            vec![
                vec![
                    ("Query_1".to_string(), b"ACGT".to_vec(), place(4, 10)),
                    ("Query_3".to_string(), b"GGGG".to_vec(), place(4, 8)),
                ],
                // `d` starts just past its end: an interval without letters; `e` keeps one letter.
                vec![
                    ("Query_4".to_string(), Vec::new(), place(4, 4)),
                    ("Query_5".to_string(), b"A".to_vec(), place(4, 5)),
                ],
            ]
        );
        assert_eq!(kept.len(), 4);
    }

    // A batch of skipped records alone is delivered empty (NCBI's query factory refuses it), also
    // at the end of the file.
    #[test]
    fn a_last_batch_of_skipped_records_is_delivered_empty() {
        let range = SequenceRange { from: 4, to: 7 };
        let (delivered, _) = batches(b">a\nACGTACGTAC\n>b\nACG\n>c\nA\n", &range, 10);
        assert_eq!(delivered.len(), 2);
        assert_eq!(delivered[0].len(), 1);
        assert!(delivered[1].is_empty());
        let (delivered, kept) = batches(b">a\nACG\n", &range, 10);
        assert!(kept.is_empty());
        assert_eq!(delivered, vec![Vec::new()]);
    }

    // A record whose letter before the interval is lower case keeps its whole letters after the
    // cut (that letter's run, which a split query's chunks restrict), whether the run ends there
    // or goes on into the interval; a record with runs only elsewhere, or none, keeps only its
    // interval.
    #[test]
    fn a_cut_record_keeps_its_whole_letters_for_a_lower_case_letter_before_it() {
        let range = SequenceRange { from: 4, to: 7 };
        let data = b">a\nACGtACGTAC\n>b\nAcgtaCGTAC\n>c\naCGTaCGTac\n>d\nACGTACGTAC\n";
        let (_, kept) = batches(data, &range, 100);
        assert_eq!(kept[0].sequence, b"ACGT");
        assert!(kept[0].lowercase_masks.is_empty());
        assert_eq!(kept[0].archive_letters(), b"ACGtACGTAC");
        assert_eq!(kept[0].archive_lowercase_runs(), vec![(3, 3)]);
        assert_eq!(kept[1].lowercase_masks, vec![(0, 0)]);
        assert_eq!(kept[1].archive_letters(), b"AcgtaCGTAC");
        assert_eq!(kept[2].lowercase_masks, vec![(0, 0)]);
        assert_eq!(kept[2].archive_letters(), b"ACGT");
        assert_eq!(kept[3].archive_letters(), b"ACGT");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1041-1046
    // ```c
    //       seqloc->ssr->left = MAX(0, seqloc->ssr->left - from);
    //       seqloc->ssr->right = MIN(seqloc->ssr->right, to) - from;
    //       /* If this mask location does not intersect the [from,to] interval,
    //          do not add it to the newly constructed list and free its contents. */
    // ```
    #[test]
    fn a_cut_clips_and_moves_the_lower_case_masks() {
        let mut record = parse_fasta(b">q\nacgTTTTaaaCCgg\n", false, true)
            .unwrap()
            .remove(0);
        assert_eq!(record.lowercase_masks, vec![(0, 2), (7, 9), (12, 13)]);
        let mut kept = record.clone();
        record.cut_in_place(2, 9, false);
        assert_eq!(record.sequence, b"GTTTTAA");
        assert_eq!(record.lowercase_masks, vec![(0, 0), (5, 6)]);
        assert_eq!(
            record.placement(),
            RecordPlacement {
                offset: 2,
                length: 14
            }
        );
        // Only a record that keeps its whole letters has them (and its own runs) after the cut.
        assert_eq!(record.archive_letters(), b"GTTTTAA");
        assert_eq!(record.archive_lowercase_runs(), vec![(0, 0), (5, 6)]);
        kept.cut_in_place(2, 9, true);
        assert_eq!(kept.sequence, record.sequence);
        assert_eq!(kept.lowercase_masks, record.lowercase_masks);
        assert_eq!(kept.placement(), record.placement());
        assert_eq!(kept.archive_letters(), b"acgTTTTaaaCCgg");
        assert_eq!(
            kept.archive_lowercase_runs(),
            vec![(0, 2), (7, 9), (12, 13)]
        );
        let whole = parse_fasta(b">q\nACGT\n", false, true).unwrap().remove(0);
        assert_eq!(
            whole.placement(),
            RecordPlacement {
                offset: 0,
                length: 4
            }
        );
    }
}

#[cfg(test)]
mod defline_tests {
    use super::*;
    use std::io::Cursor;

    /// The records of `data` read as a nucleotide query file, with or without `-parse_deflines`.
    fn read(data: &[u8], parse_deflines: bool) -> Result<Vec<FastaRecord>> {
        let mut stream = FastaStream::new(Cursor::new(data.to_vec()));
        stream.parse_deflines = parse_deflines;
        parse_fasta_stream_keeping(
            &mut stream,
            false,
            false,
            100,
            false,
            None,
            &mut |_| Ok(()),
            &mut |_| Ok(()),
        )
    }

    /// `read` with NCBI's data loaders on (`FastaStream::seq_id_role`): the titles and letters.
    fn read_with_loaders(data: &[u8], parse_deflines: bool) -> Result<Vec<(String, Vec<u8>)>> {
        let mut stream = FastaStream::new(Cursor::new(data.to_vec()));
        stream.parse_deflines = parse_deflines;
        stream.seq_id_role = Some("query");
        let records = parse_fasta_stream_keeping(
            &mut stream,
            false,
            false,
            100,
            false,
            None,
            &mut |_| Ok(()),
            &mut |_| Ok(()),
        )?;
        Ok(records
            .into_iter()
            .map(|record| (record.title, record.sequence))
            .collect())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:130-133
    // ```c++
    //     virtual CRef<CSeq_entry> ReadOneSeq(ILineErrorListener * pMessageListener) {
    //
    //         const string line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
    //         if ( !line.empty() && isalnum(line.data()[0]&0xff) ) {
    // ```
    // With the data loaders on, BLASTX's reader refuses the first lines that the other programs'
    // reader refuses, with that reader's text, also under `-parse_deflines` (the Seq-id is tried
    // before `CFastaReader` reads); a first line that NCBI reads as FASTA, and a Seq-id after the
    // first line, are read as with the data loaders off.
    #[test]
    fn a_first_line_ncbi_tries_as_a_seq_id_is_refused_as_the_shared_reader_refuses_it() {
        use crate::blastinput::fasta_reader::{FastaInputSource, ReadError, ReaderConfig};
        // The first lines, each before a record with LF line ends, and whole inputs: CR-only and
        // CRLF line ends, and a file of the line alone without a line end.
        let first_lines: [&[u8]; 9] = [
            b"lcl|abc",
            b"12345",
            b"AB123456",
            b"P12345",
            b"gnl|db|x",
            b"  lcl|abc  ",
            b"lcl|abc\r",
            b"\x0blcl|abc",
            b"\x0c12345",
        ];
        let inputs = first_lines
            .iter()
            .map(|first| [*first, b"\n>q1 x\nACGTACGT\n"].concat())
            .chain([
                b"lcl|abc\r>q1 x\rACGTACGT\r".to_vec(),
                b"12345\r\n>q1 x\r\nACGTACGT\r\n".to_vec(),
                b"lcl|abc".to_vec(),
            ]);
        for data in inputs {
            let config = ReaderConfig::query("BLASTX", false, true);
            let expected =
                match FastaInputSource::from_bytes(&data, config).next_sequence(&mut |_| Ok(())) {
                    Err(ReadError::Unsupported(refusal)) => refusal.to_string(),
                    other => panic!("{data:?}: {other:?}"),
                };
            assert!(
                expected.contains("not supported by LOSAT's BLASTX"),
                "{expected}"
            );
            for parse_deflines in [false, true] {
                let error = read_with_loaders(&data, parse_deflines).unwrap_err();
                assert_eq!(error.to_string(), expected, "{data:?}");
            }
            // With the data loaders off the line is data: the letters of a first record without a
            // defline, or `CheckDataLine`'s error for a line that is not plausible data (`12345`).
            match read(&data, false) {
                Ok(records) => assert!(records[0].title.is_empty(), "{data:?}"),
                Err(error) => assert!(error.is::<FastaParseError>(), "{data:?}: {error}"),
            }
        }
        let fasta: [&[u8]; 5] = [
            b"ACGTACGTAC\n>q1 x\nACGT\n",
            b"ACGT ACGT\n>q1 x\nACGT\n",
            b"\nlcl|abc\n>q1 x\nACGT\n",
            b";c\nlcl|abc\n>q1 x\nACGT\n",
            b">lcl|abc x\nACGT\n>q1 x\nACGT\n",
        ];
        for data in fasta {
            let without: Vec<_> = read(data, false)
                .unwrap()
                .into_iter()
                .map(|record| (record.title, record.sequence))
                .collect();
            assert_eq!(read_with_loaders(data, false).unwrap(), without, "{data:?}");
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:383-398
    // ```c++
    //         if (c == '!'  ||  c == '#' || c == ';') {
    //             // no content, just a comment or blank line
    //             continue;
    //         } else if (need_defline) {
    //             if (TestFlag(fDLOptional)) {
    //                 ParseDefLine(">", pMessageListener);
    //                 need_defline = false;
    //             } else {
    //                 const auto lineNum = LineNumber();
    //                 GetLineReader().UngetLine();
    //                 NCBI_THROW2(CObjReaderParseException, eNoDefline,
    //                             "CFastaReader: Input doesn't start with"
    //                             " a defline or comment around line " + NStr::NumericToString(lineNum),
    //                              lineNum);
    //             }
    //         }
    // ```
    // Under -parse_deflines a gap line before any defline stops the read at its line, whatever
    // its size or modifiers; without it, the gap line starts a record without a defline. A gap
    // line inside a record is that record's letters either way.
    #[test]
    fn a_gap_line_before_any_defline_stops_a_parse_deflines_read() {
        for (data, line) in [
            (&b">?5\nACGT\n>r1\nACGT\n"[..], 1),
            (b">?5", 1),
            (b">?0\nACGT\n", 1),
            (b">?abc\nACGT\n", 1),
            (b">?unk100\n>?7\nACGT\n", 1),
            (b"  >?5\r\nACGT\r\n", 1),
            (b";c\n\n#x\n>?5\nACGT\n", 4),
        ] {
            let error = read(data, true).unwrap_err();
            assert_eq!(
                error.downcast_ref::<NoDeflineError>().map(|e| (e.0, e.to_string())),
                Some((
                    line,
                    format!(
                        "CFastaReader: Input doesn't start with a defline or comment around line {line}"
                    )
                )),
                "{}",
                String::from_utf8_lossy(data)
            );
        }
        let records = read(b">?5\nACGT\n>r1\nACGT\n", false).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].sequence, b"NNNNNACGT");
        let records = read(b">r0\nAC\n>?5\nGT\n", true).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, b"ACNNNNNGT");
    }
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:774-780
    // ```c++
    // void CFastaReader::ParseDataLine(
    //     const TStr& s, ILineErrorListener * pMessageListener)
    // {
    //     if( NStr::StartsWith(s, ">?") ) {
    //         ParseGapLine(s, pMessageListener);
    //         return;
    //     }
    // ```
    // BLASTX's reader reads a gap line as the reader of the other programs reads it
    // (`fasta_reader::reader::gap_line_size`): the same letters, lower-case runs and warnings,
    // whatever the size and the modifiers, before any defline, inside a record, as its last line
    // and inside a lower-case run. NCBI 2.17.0 reads them so in all five programs (probe
    // `ncbi-parity-20261019/SM/gapline`).
    #[test]
    fn a_gap_line_reads_as_the_reader_of_the_other_programs_reads_it() {
        use crate::blastinput::fasta_reader::{
            FastaInputSource, ParseErrorCode, ReadError, ReaderConfig,
        };
        let gaps: [&[u8]; 23] = [
            b">?5",
            b">? 5 ",
            b">?0",
            b">?abc",
            b">?-5",
            b">?+5",
            b">?5abc",
            b">?",
            b">?unk",
            b">?unk7",
            b">?unk 7",
            b">?99999999999",
            b">?4294967296",
            b">?5 [gap-type=short-arm]",
            b">?5 [gap-type=bogus]",
            b">?5 [linkage-evidence=paired-ends]",
            b">?5 [gap-type=within scaffold][linkage-evidence=paired-ends]",
            b">?5 [gap-type=short-arm][gap-type=centromere]",
            b">?5 [foo=bar]",
            b">?5 [gap-type",
            b">?0 [gap-type=short-arm] x",
            b">?5 [gap-type=caf\xe9]",
            b">?5 [gap-type=a\0b] [linkage-evidence=pcr]",
        ];
        // BLASTX's records of `input` (a nucleotide query file, `-lcase_masking`): the letters with
        // the lower-case runs in lower case, as the reader of the other programs gives them, and
        // the warnings as BLASTX writes them.
        let blastx = |input: &[u8]| {
            let mut warnings = Vec::new();
            let records = parse_fasta_stream_keeping(
                &mut FastaStream::new(Cursor::new(input.to_vec())),
                false,
                true,
                100,
                false,
                None,
                &mut |_| Ok(()),
                &mut |warning| {
                    crate::algorithm::blastx::native::reader_warning(warning, &mut warnings)
                },
            )
            .unwrap();
            let letters: Vec<Vec<u8>> = records
                .iter()
                .map(|record| {
                    let mut letters = record.sequence.to_ascii_uppercase();
                    for &(left, right) in &record.lowercase_masks {
                        letters[left as usize..=right as usize].make_ascii_lowercase();
                    }
                    letters
                })
                .collect();
            (letters, warnings)
        };
        for gap in gaps {
            for input in [
                [gap, b"\nACGT\n>r1\nAC\n"].concat(),
                [b">r0\nAC\n", gap, b"\nGT\n>r1\nAC\n"].concat(),
                [b">r0\nAC\n", gap, b"\n>r1\nAC\n"].concat(),
                [b">r0\nac\n", gap, b"\ngtAC\n>r1\nac\n>?2\n"].concat(),
            ] {
                let context = String::from_utf8_lossy(&input).into_owned();
                let (letters, warnings) = blastx(&input);
                let mut source = FastaInputSource::from_bytes(
                    &input,
                    ReaderConfig::query("BLASTX", false, false),
                );
                let (mut messages, mut expected) = (Vec::new(), Vec::new());
                while !source.end() {
                    match source.next_sequence(&mut |message: &[u8]| {
                        messages.extend_from_slice(message);
                        Ok(())
                    }) {
                        Ok(record) => expected.push(record.sequence),
                        Err(ReadError::Parse {
                            code: ParseErrorCode::Eof,
                            ..
                        }) => break,
                        Err(error) => panic!("{context}: {error}"),
                    }
                }
                assert_eq!(letters, expected, "{context}");
                assert_eq!(warnings, messages, "{context}");
            }
        }
        // NCBI 2.17.0's warnings and letters (`blastx -query` with `>?abc` inside a record, probe
        // case `blastx.q.mid_abc.f6.sub`; `>?5 [gap-type=caf\xe9]` keeps its byte).
        let (letters, warnings) = blastx(b">r0\nAC\n>?abc\nGT\n>r1\nAC\n>?5 [gap-type=caf\xe9]\n");
        assert_eq!(letters, vec![b"ACNGT".to_vec(), b"ACNNNNN".to_vec()]);
        assert_eq!(
            warnings,
            b"CFastaReader: Bad gap size at line 3\nCFastaReader: Problem parsing gap mods at line 3\nUnknown gap-type: caf\xe9\n"
        );
    }
}

#[cfg(test)]
mod query_settings_tests {
    use super::*;

    /// The directory of this test's files (removed at its end).
    fn dir() -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "losat-blastx-query-settings-{}",
            std::process::id()
        ))
    }

    /// A file of this test with the bytes `text`, open for reading.
    fn file(tag: &str, text: &[u8]) -> std::fs::File {
        std::fs::create_dir_all(dir()).unwrap();
        let path = dir().join(tag);
        std::fs::write(&path, text).unwrap();
        std::fs::File::open(&path).unwrap()
    }

    /// The records of `stream`, of the query (`protein` false) or the subject role.
    fn read(stream: &mut FastaStream<std::fs::File>, protein: bool) -> Result<Vec<FastaRecord>> {
        parse_fasta_stream_keeping(
            stream,
            protein,
            false,
            100,
            false,
            None,
            &mut |_| Ok(()),
            &mut |_| Ok(()),
        )
    }

    // The read before the search reads as the search does (under `-parse_deflines` a first line
    // that is not a defline is the reader's error) and adds its ids to a copy of the search's: a
    // subject with the id of a query that only that read has read is not refused, and that read
    // stops at a query with the id of a subject the search has read, as the search does.
    #[test]
    fn a_read_before_the_search_reads_as_the_search_and_keeps_the_search_ids() {
        let settings = QueryStreamSettings {
            parse_deflines: true,
            id_scope: Some(Default::default()),
            strand: crate::report::archive::bioseq::SeqLocStrand::Both,
            data_loaders: false,
        };
        let error = read(
            &mut settings.stream_before_search(file("nodef", b"ACGT\n>q1 x\nACGT\n")),
            false,
        )
        .unwrap_err();
        assert!(error.is::<NoDeflineError>(), "{error}");
        let queries = read(
            &mut settings.stream_before_search(file("query", b">q1 x\nACGTACGT\n")),
            false,
        )
        .unwrap();
        assert_eq!(queries.len(), 1);
        read(
            &mut settings.stream(file("subject", b">q1 y\nMKVL\n")),
            true,
        )
        .unwrap();
        let error = read(
            &mut settings.stream_before_search(file("again", b">q1 x\nACGTACGT\n")),
            false,
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("a query and a subject record have the Seq-id lcl|q1"),
            "{error}"
        );
        std::fs::remove_dir_all(dir()).unwrap();
    }
}
