//! NCBI BLAST+'s reading of FASTA input for BLASTN, TBLASTX, TBLASTN and BLASTP: the
//! input source (`CBlastFastaInputSource`, `CBlastInput`), the readers
//! (`CBlastInputReader`, `CCustomizedFastaReader`, `CFastaReader`) and the line reader
//! (`CStreamLineReader`). Nucleotide input uses `fAssumeNuc`, protein input `fAssumeProt`;
//! the other flags are the same.
//!
//! A record keeps what the search and the reports use: the generated local ID
//! (`Query_N`, `Subject_N`), the title bytes, the residues, and the reader's messages. The
//! lowercase mask (`x_OpenMask`/`x_CloseMask`) is kept as the case of the residues: a
//! residue in a masked interval is in lower case, every other residue in upper case.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:316-361
//! ```c++
//! void
//! CBlastFastaInputSource::x_InitInputReader()
//! {
//!     CFastaReader::TFlags flags = m_Config.GetBelieveDeflines() ?
//!                                     CFastaReader::fParseRawID:
//!                                     (CFastaReader::fNoParseID |
//!                                      CFastaReader::fDLOptional);
//!
//!     // Allow CFastaReader fSkipCheck flag to be set based
//!     // on new CBlastInputSourceConfig property - GetSkipSeqCheck() -RMH-
//!     flags += ( m_Config.GetSkipSeqCheck() ? CFastaReader::fSkipCheck : 0 );
//!
//!     flags += (m_ReadProteins
//!               ? CFastaReader::fAssumeProt
//!               : CFastaReader::fAssumeNuc);
//!     const char* env_var = getenv("BLASTINPUT_GEN_DELTA_SEQ");
//!     if (env_var == NULL || (env_var && string(env_var) == kEmptyStr)) {
//!         flags += CFastaReader::fNoSplit;
//!     }
//!     // This is necessary to enable the ignoring of gaps in classes derived from
//!     // CFastaReader
//!
//!    	flags+= CFastaReader::fHyphensIgnoreAndWarn;
//!
//!     flags+= CFastaReader::fDisableNoResidues;
//!     // Do not check more than few characters in local ID for illegal characters.
//!     // Illegal characters can be things like = and we want to let those through.
//!     flags+= CFastaReader::fQuickIDCheck;
//!
//!     if (m_Config.GetDataLoaderConfig().UseDataLoaders()) {
//!         m_InputReader.reset
//!             (new CBlastInputReader(m_Config.GetDataLoaderConfig(),
//!                                    m_ReadProteins,
//!                                    m_Config.RetrieveSeqData(),
//!                                    m_Config.GetSeqLenThreshold2Guess(),
//!                                    *m_LineReader,
//!                                    flags));
//!     } else {
//!         m_InputReader.reset(new CCustomizedFastaReader(*m_LineReader, flags,
//!                                        m_Config.GetSeqLenThreshold2Guess()));
//!     }
//! ```
//! `-parse_deflines` (`fParseRawID`) is rejected by LOSAT's programs; `fSkipCheck` is not
//! set by the BLAST programs.

mod reader;
mod seq_id;
mod stream;
#[cfg(test)]
mod tests;

use std::borrow::Cow;
use std::io::{Read, Seek};

pub(crate) use stream::FastaStream;
use stream::LineReader;

/// How an input is read: the molecule, the role's local ID prefix and whether the data
/// loaders are configured (the Seq-id path of `CBlastInputReader`).
#[derive(Clone, Copy, Debug)]
pub struct ReaderConfig {
    /// `fAssumeProt` (protein input) or `fAssumeNuc`.
    pub protein: bool,
    /// `SetQueryLocalIdMode` (`Query_`) or `SetSubjectLocalIdMode` (`Subject_`).
    pub id_prefix: &'static str,
    /// `SDataLoaderConfig::UseDataLoaders()`: `CBlastInputReader` tries the first line of
    /// a record as a Seq-id.
    pub data_loaders: bool,
    /// The program and role named in LOSAT's rejections.
    pub program: &'static str,
    pub role: &'static str,
}

impl ReaderConfig {
    /// The query reader of `program`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:66-74
    /// ```c++
    ///     // Set an appropriate default for the strand
    ///     if (m_Strand == eNa_strand_other) {
    ///         m_Strand = (m_DLConfig.m_IsLoadingProteins)
    ///             ? eNa_strand_unknown
    ///             : eNa_strand_both;
    ///     }
    ///     SetQueryLocalIdMode();
    /// }
    /// ```
    pub fn query(program: &'static str, protein: bool, data_loaders: bool) -> Self {
        Self {
            protein,
            id_prefix: "Query_",
            data_loaders,
            program,
            role: "query",
        }
    }

    /// The subject reader of `program`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:230-237
    /// ```c++
    ///     SDataLoaderConfig dlconfig(read_proteins);
    ///     dlconfig.OptimizeForWholeLargeSequenceRetrieval();
    ///
    ///     CBlastInputSourceConfig iconfig(dlconfig);
    ///     iconfig.SetRange(range);
    ///     iconfig.SetBelieveDeflines(parse_deflines);
    ///     iconfig.SetLowercaseMask(use_lcase_masking);
    ///     iconfig.SetSubjectLocalIdMode();
    /// ```
    pub fn subject(program: &'static str, protein: bool, data_loaders: bool) -> Self {
        Self {
            protein,
            id_prefix: "Subject_",
            data_loaders,
            program,
            role: "subject",
        }
    }
}

/// One record as NCBI's FASTA input source reads it.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct FastaRecord {
    /// The generated local ID (`lcl|Query_N`, `lcl|Subject_N`, without `lcl|`).
    pub local_id: String,
    /// The title: the defline after `>` and the white space that follows it, up to the
    /// first byte below 0x20, without white space at its ends (`x_ApplyMods`); empty when
    /// the record has none. Any other bytes, including bytes that are not UTF-8.
    pub title: Vec<u8>,
    /// The residues (IUPAC letters; nucleotide `U` as `T`; `>?` gaps as `N` or `X`), in
    /// lower case where the lowercase mask covers them.
    pub sequence: Vec<u8>,
    /// The reader's messages written while the record was read, in order (each one a
    /// line; NCBI writes them to standard error as it reads).
    pub warnings: Vec<u8>,
}

impl FastaRecord {
    /// A record from its parts (tests and the ABI's records).
    pub fn new(local_id: impl Into<String>, title: &[u8], sequence: &[u8]) -> Self {
        Self {
            local_id: local_id.into(),
            title: title.to_vec(),
            sequence: sequence.to_vec(),
            warnings: Vec::new(),
        }
    }

    /// The residues.
    pub fn seq(&self) -> &[u8] {
        &self.sequence
    }

    /// The title bytes (empty without one).
    pub fn title(&self) -> &[u8] {
        &self.title
    }

    /// The ID that the tabular formats and the reports show for a record with a local
    /// ID: the title up to its first space, or the local ID when there is no title.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:474-504
    /// (`s_ReplaceLocalId`; see `report::` for its callers)
    pub fn shown_id(&self) -> &[u8] {
        if self.title.is_empty() {
            self.local_id.as_bytes()
        } else {
            let end = self
                .title
                .iter()
                .position(|&byte| byte == b' ')
                .unwrap_or(self.title.len());
            &self.title[..end]
        }
    }

    /// The record cut to the residues `from..to_exclusive` of a range (`-query_loc`,
    /// `-subject_loc`): the same local ID and title, the letters (their case too) of the
    /// interval, and no messages (the reader wrote them when it read the whole record).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:458-466
    /// ```c++
    ///     // set sequence range
    ///     retval->SetInt().SetFrom(from);
    ///     retval->SetInt().SetTo((to > 0 && to < seqlen) ? to : (seqlen-1));
    ///
    ///     // set ID
    ///     retval->SetInt().SetId().Assign(*FindBestChoice(itr->GetId(), CSeq_id::BestRank));
    ///
    ///     return retval;
    /// }
    /// ```
    /// NCBI keeps the whole record and searches the interval of its location; LOSAT keeps
    /// the interval's letters and where it lies (`seq_range::Placements`).
    pub fn cut(&self, from: usize, to_exclusive: usize) -> Self {
        Self {
            local_id: self.local_id.clone(),
            title: self.title.clone(),
            sequence: self.sequence[from..to_exclusive].to_vec(),
            warnings: Vec::new(),
        }
    }

    /// The record of a `bio` record (the bridge of the port plan, §1): only where LOSAT's
    /// checks guarantee that `bio` reads the input as NCBI's reader does (ABI v1, plan TD-1;
    /// the Web adapter until step S10). `n` is the `N` of the local ID (`Query_N`,
    /// `Subject_N`, 1-based) and `prefix` its prefix; `protein` is the molecule of the input
    /// (`fAssumeProt`, else `fAssumeNuc`).
    ///
    /// The title is the defline as `bio` splits it (the ID, a space and the description),
    /// without the white space at its start, which NCBI's defline parser skips (`bio` gives
    /// such a defline an empty ID); that is NCBI's title for the deflines those checks
    /// accept. The messages are the title warning that NCBI writes at the end of such a
    /// record.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:209-213
    /// ```c++
    ///     // trim leading whitespace from title (is this appropriate?)
    ///     while (title_start < len
    ///         &&  isspace((unsigned char)defline[title_start])) {
    ///         ++title_start;
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2038-2043
    /// ```c++
    ///     NStr::TruncateSpacesInPlace(processed_title);
    ///     if (!processed_title.empty()) {
    ///         auto pDesc = Ref(new CSeqdesc());
    ///         pDesc->SetTitle() = processed_title;
    ///         bioseq.SetDescr().Set().push_back(std::move(pDesc));
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:469-482
    /// ```c++
    ///     auto n = m_Counter.load();
    ///     if (advance)
    ///         m_Counter++;
    ///
    ///     if (m_Prefix.empty()  &&  m_Suffix.empty()) {
    ///         seq_id->SetLocal().SetId(n);
    ///     } else {
    ///         string& id = seq_id->SetLocal().SetStr();
    ///         id.reserve(128);
    ///         id += m_Prefix;
    ///         id += NStr::IntToString(n);
    ///         id += m_Suffix;
    ///     }
    ///     return seq_id;
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1444-1447
    /// ```c++
    ///         CRef<CSeq_data> data(new CSeq_data(m_SeqData, format));
    ///         if ( !TestFlag(fLeaveAsText) ) {
    ///             CSeqportUtil::Pack(data, inst.GetLength());
    ///         }
    /// ```
    /// A nucleotide `U` is stored as `T` (`u` as `t` inside the lowercase mask), as the
    /// reader stores it (`reader.rs`, `assemble_seq`).
    pub fn from_bio(
        record: &bio::io::fasta::Record,
        n: usize,
        prefix: &str,
        protein: bool,
    ) -> Self {
        let mut title = Vec::with_capacity(record.id().len() + 1);
        title.extend_from_slice(record.id().as_bytes());
        if let Some(desc) = record.desc() {
            title.push(b' ');
            title.extend_from_slice(desc.as_bytes());
        }
        // `isspace` of the C locale.
        let start = title
            .iter()
            .position(|&byte| !matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'))
            .unwrap_or(title.len());
        title.drain(..start);
        let mut sequence = record.seq().to_vec();
        if !protein {
            for residue in sequence.iter_mut() {
                *residue = match *residue {
                    b'U' => b'T',
                    b'u' => b't',
                    other => other,
                };
            }
        }
        let mut warnings = Vec::new();
        if let Some(warning) = reader::seq_data_in_title_warning(&title, protein) {
            warnings.extend_from_slice(warning);
            warnings.push(b'\n');
        }
        Self {
            local_id: format!("{prefix}{n}"),
            title,
            sequence,
            warnings,
        }
    }
}

/// A record as the programs use it: its residues, the record cut to a range, and its title
/// bytes (the ranges in `seq_range.rs`, the query warnings, the masks and the lookup
/// tables). The bridge of the port plan (§1) also implemented it for the `bio` records until
/// every program read with NCBI's reader (step S8).
pub trait InputRecord: Sized {
    /// The residues (the search's letters, lower case where the lowercase mask covers
    /// them).
    fn seq(&self) -> &[u8];

    /// The record cut to the residues `from..to_exclusive` (`FastaRecord::cut`).
    fn cut(&self, from: usize, to_exclusive: usize) -> Self;

    /// The title (`CBlastQuerySourceOM::GetTitle`: the record's title descriptor; empty
    /// without one).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:366-372
    /// ```c++
    ///     string title(kEmptyStr);
    ///     if (bh.CanGetDescr())
    ///     {
    ///     	const CSeq_descr::Tdata& descr = bh.GetDescr();
    ///     	ITERATE(CSeq_descr::Tdata, desc, descr) {
    ///         	if ((*desc)->Which() == CSeqdesc::e_Title && title == kEmptyStr) {
    ///             		title = (*desc)->GetTitle();
    /// ```
    fn title_bytes(&self) -> Cow<'_, [u8]>;
}

impl InputRecord for FastaRecord {
    fn seq(&self) -> &[u8] {
        &self.sequence
    }

    fn cut(&self, from: usize, to_exclusive: usize) -> Self {
        FastaRecord::cut(self, from, to_exclusive)
    }

    fn title_bytes(&self) -> Cow<'_, [u8]> {
        Cow::Borrowed(&self.title)
    }
}

/// `CObjReaderParseException`'s codes that the reader throws with these flags.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ParseErrorCode {
    /// `eEOF`: the batch ends (`GetNextSeqBatch`, `GetAllSeqs`).
    Eof,
    /// `eFormat` (`CheckDataLine`).
    Format,
}

/// Why a record was not read.
#[derive(Debug)]
pub enum ReadError {
    /// A `CObjReaderParseException` with its message (`GetMsg()`) and line.
    Parse {
        code: ParseErrorCode,
        message: String,
        line: u64,
    },
    /// Input that LOSAT rejects explicitly (`not supported by LOSAT's <PROGRAM>`).
    Unsupported(anyhow::Error),
    /// The warnings could not be written.
    Write(std::io::Error),
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadError::Parse { message, .. } => f.write_str(message),
            ReadError::Unsupported(error) => write!(f, "{error:#}"),
            ReadError::Write(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for ReadError {}

impl ReadError {
    /// The error as the programs end with it: a reader exception is NCBI's `BLAST query
    /// error: <message>` with exit 1, a rejection is LOSAT's, and a failed write of the
    /// messages is the write error.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:181-184
    /// ```c++
    ///     catch (const CObjReaderParseException& e) {                             \
    ///         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
    ///         exit_code = BLAST_INPUT_ERROR;                                      \
    ///     }                                                                       \
    /// ```
    pub fn into_app_error(self) -> anyhow::Error {
        match self {
            ReadError::Parse { message, .. } => crate::cli::NativeError {
                exit: 1,
                message: format!("BLAST query error: {message}\n"),
            }
            .into(),
            ReadError::Unsupported(error) => error,
            ReadError::Write(error) => error.into(),
        }
    }
}

/// NCBI's `CBlastFastaInputSource` over one opened input.
pub struct FastaInputSource<R: Read> {
    reader: reader::FastaReader<R>,
}

impl FastaInputSource<std::fs::File> {
    /// An input source over an opened file or standard input (`CStreamLineReader` over
    /// the argument's stream).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:291-300
    /// ```c++
    /// CBlastFastaInputSource::CBlastFastaInputSource(CNcbiIstream& infile,
    ///                                        const CBlastInputSourceConfig& iconfig)
    ///     : m_Config(iconfig),
    ///       m_LineReader(iconfig.GetConvertGapsToNs() ?
    ///                    new CStreamLineReaderConverter(infile) :
    ///                    new CStreamLineReader(infile)),
    ///       m_ReadProteins(iconfig.IsProteinInput())
    /// {
    ///     x_InitInputReader();
    /// }
    /// ```
    /// `GetConvertGapsToNs` is set only by the mapper (magicblast).
    pub fn from_file(file: std::fs::File, config: ReaderConfig) -> Self {
        Self::from_stream(FastaStream::from_file(file), config)
    }

    /// An input source over standard input (`-`, NCBI's `cin`; `FastaStream::from_standard_input`).
    pub fn from_standard_input(file: std::fs::File, config: ReaderConfig) -> Self {
        Self::from_stream(FastaStream::from_standard_input(file), config)
    }

    /// The input source of an opened `-query` or `-subject` argument: standard input for `-`
    /// (`from_standard_input`), else a file (`from_file`).
    pub fn from_argument(
        path: &std::path::Path,
        file: std::fs::File,
        config: ReaderConfig,
    ) -> Self {
        if path.as_os_str() == "-" {
            Self::from_standard_input(file, config)
        } else {
            Self::from_file(file, config)
        }
    }
}

impl<'a> FastaInputSource<&'a [u8]> {
    /// An input source over bytes in memory (the ABI's registered inputs), read as a
    /// regular file with the same bytes (`FastaStream::from_bytes`; test
    /// `bytes_read_as_a_file_with_the_same_bytes`).
    pub fn from_bytes(bytes: &'a [u8], config: ReaderConfig) -> Self {
        Self::from_stream(FastaStream::from_bytes(bytes), config)
    }
}

impl<'a> FastaInputSource<&'a [u8]> {
    /// `from_bytes` with LOSAT's bulk paths off: one byte at a time, the reference that
    /// the tests compare the bulk paths with.
    #[cfg(test)]
    pub(crate) fn from_bytes_byte_at_a_time(bytes: &'a [u8], config: ReaderConfig) -> Self {
        let mut stream = FastaStream::from_bytes(bytes);
        stream.bulk = false;
        Self::from_stream(stream, config)
    }
}

impl<R: Read + Seek> FastaInputSource<R> {
    /// `IsIStreamEmpty` on the input, which the programs call on the query's stream before
    /// they build the input source: true when the input has only white space (`Query is
    /// Empty!`), false for a stream without a position, such as a pipe, whatever it holds.
    /// Call it before the first record is read.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:845-875
    /// ```c++
    /// bool
    /// IsIStreamEmpty(CNcbiIstream & in)
    /// {
    /// #ifdef NCBI_OS_MSWIN
    /// 	char c;
    /// 	in.setf(ios::skipws);
    /// 	if (!(in >> c))
    /// 		return true;
    /// 	in.unget();
    /// 	return false;
    /// #else
    /// 	char c;
    /// 	CNcbiStreampos orig_p = in.tellg();
    /// 	// Piped input
    /// 	if(orig_p < 0)
    /// 		return false;
    ///
    /// 	IOS_BASE::iostate orig_state = in.rdstate();
    /// 	IOS_BASE::fmtflags orig_flags = in.setf(ios::skipws);
    ///
    /// 	if(! (in >> c))
    /// 		return true;
    ///
    /// 	in.seekg(orig_p);
    /// 	in.flags(orig_flags);
    /// 	in.clear();
    /// 	in.setstate(orig_state);
    ///
    /// 	return false;
    /// #endif
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:209-214
    /// ```c++
    ///         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())) {
    ///            	ERR_POST(Warning << "Query is Empty!");
    ///            	return BLAST_EXIT_SUCCESS;
    ///         }
    ///         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
    ///         CBlastInput input(&fasta);
    /// ```
    /// LOSAT takes the Linux branch on every platform: the oracle is NCBI BLAST+ on Linux
    /// (`AUTHORITY.md` §C1).
    pub fn stream_is_empty(&mut self) -> bool {
        self.reader.lines.stream_is_empty()
    }
}

impl<R: Read> FastaInputSource<R> {
    fn from_stream(stream: FastaStream<R>, config: ReaderConfig) -> Self {
        Self {
            reader: reader::FastaReader::new(LineReader::new(stream), config),
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:370-374
    /// ```c++
    /// bool
    /// CBlastFastaInputSource::End()
    /// {
    ///     return m_LineReader->AtEOF();
    /// }
    /// ```
    pub fn end(&mut self) -> bool {
        self.reader.lines.at_eof()
    }

    /// The next record (`GetNextSequence`, without the range and strand of
    /// `x_FastaToSeqLoc`, which the programs apply: `seq_range.rs`). The reader's messages
    /// go to `warn` as NCBI writes them, and are also kept in the record.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:376-387
    /// ```c++
    /// CRef<CSeq_loc>
    /// CBlastFastaInputSource::x_FastaToSeqLoc(CRef<objects::CSeq_loc>& lcase_mask,
    ///                                         CScope& scope)
    /// {
    ///     static const TSeqRange kEmptyRange(TSeqRange::GetEmpty());
    ///     CRef<CBlastScopeSource> query_scope_source;
    ///
    ///     if (m_Config.GetLowercaseMask())
    ///         lcase_mask = m_InputReader->SaveMask();
    ///
    ///     CRef<CSeq_entry> seq_entry(m_InputReader->ReadOneSeq());
    /// ```
    /// The molecule checks after it (`Nucleotide FASTA provided for protein sequence` and
    /// the reverse) cannot fail: `AssignMolType` sets the molecule of the input.
    ///
    /// A line that NCBI's data loaders would fetch as a Seq-id gives
    /// `ReadError::Unsupported`; the source is then past that one line, with the local-ID
    /// counter unchanged, so a further call reads from the next line, as NCBI's next
    /// `ReadOneSeq` does. The programs stop at the rejection (`read_all`, `read_queries`),
    /// where NCBI skips a query it cannot fetch and fails on such a subject.
    pub fn next_sequence(
        &mut self,
        warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
    ) -> Result<FastaRecord, ReadError> {
        let mut warnings = Vec::new();
        let mut record = self.reader.read_one_seq(&mut |message: &[u8]| {
            warnings.extend_from_slice(message);
            warn(message)
        })?;
        record.warnings = warnings;
        Ok(record)
    }
}

/// `IsIStreamEmpty` on an opened input (true when it has only white space; false for a
/// stream without a position, such as a pipe).
pub fn input_is_empty<R: Read + Seek>(stream: &mut FastaStream<R>) -> bool {
    stream::stream_is_empty(stream)
}

/// All records of an input (`CBlastInput::GetAllSeqs`, the subjects): the messages go to
/// `warn` as they are written; a reader error stops the reading.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:198-219
/// ```c++
/// CRef<CBlastQueryVector>
/// CBlastInput::GetAllSeqs(CScope& scope)
/// {
///     CRef<CBlastQueryVector> retval(new CBlastQueryVector);
///
///     while (!End()) {
///         try { retval->AddQuery(m_Source->GetNextSequence(scope)); }
///         catch (const CObjReaderParseException& e) {
///             auto err = e.GetErrCode();
///             if (err == CObjReaderParseException::eEOF) {
///                 break;
///             } else if (err == CObjReaderParseException::eNoDefline) {
///                 CNcbiStrstream ss;
///                 ss << "Query input doesn't start with "
///                     "a defline or comment, line " << e.GetPos() << ends;
///                 NCBI_THROW(CInputException, eInvalidInput, ss.str());
///             }
///             throw;
///         }
///     }
///
///     return retval;
/// }
/// ```
/// `eNoDefline` needs a reader without `fDLOptional`.
pub fn read_all<R: Read>(
    source: &mut FastaInputSource<R>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
) -> Result<Vec<FastaRecord>, ReadError> {
    let mut records = Vec::new();
    while !source.end() {
        match source.next_sequence(warn) {
            Ok(record) => records.push(record),
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => break,
            Err(error) => return Err(error),
        }
    }
    Ok(records)
}

/// The subjects of `-subject` as NCBI reads them while it processes the arguments
/// (`ReadSequencesToBlast`): every record, with the reader's messages written to `warn` as
/// each record is read; a record whose interval of `range` (`-subject_loc`) starts more than
/// one letter past its end stops the reading after that record, with NCBI's range error. The
/// records are returned whole (`seq_range::cut_subjects` cuts them to `range`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:221-247
/// ```c++
/// CRef<CScope>
/// ReadSequencesToBlast(CNcbiIstream& in,
///                      bool read_proteins,
///                      const TSeqRange& range,
///                      bool parse_deflines,
///                      bool use_lcase_masking,
///                      CRef<CBlastQueryVector>& sequences,
///                      bool gaps_to_Ns /* = false */)
/// {
///     SDataLoaderConfig dlconfig(read_proteins);
///     dlconfig.OptimizeForWholeLargeSequenceRetrieval();
///
///     CBlastInputSourceConfig iconfig(dlconfig);
///     iconfig.SetRange(range);
///     iconfig.SetBelieveDeflines(parse_deflines);
///     iconfig.SetLowercaseMask(use_lcase_masking);
///     iconfig.SetSubjectLocalIdMode();
///     if (!read_proteins && gaps_to_Ns) {
///         iconfig.SetConvertGapsToNs(true);
///     }
///
///     CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
///     CRef<CBlastInput> input(new CBlastInput(fasta));
///     CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
///     sequences = input->GetAllSeqs(*scope);
///     return scope;
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:198-219
/// ```c++
/// CRef<CBlastQueryVector>
/// CBlastInput::GetAllSeqs(CScope& scope)
/// {
///     CRef<CBlastQueryVector> retval(new CBlastQueryVector);
///
///     while (!End()) {
///         try { retval->AddQuery(m_Source->GetNextSequence(scope)); }
///         catch (const CObjReaderParseException& e) {
///             auto err = e.GetErrCode();
///             if (err == CObjReaderParseException::eEOF) {
///                 break;
///             } else if (err == CObjReaderParseException::eNoDefline) {
///                 CNcbiStrstream ss;
///                 ss << "Query input doesn't start with "
///                     "a defline or comment, line " << e.GetPos() << ends;
///                 NCBI_THROW(CInputException, eInvalidInput, ss.str());
///             }
///             throw;
///         }
///     }
///
///     return retval;
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:450-453
/// ```c++
///     if (from > seqlen) {
///         NCBI_THROW(CInputException, eInvalidRange,
///                    "Invalid from coordinate (greater than sequence length)");
///     }
/// ```
/// `eNoDefline` needs a reader without `fDLOptional`. The range's `CInputException` is
/// `BLAST query/options error: ...` (`seq_range::record_interval`, `app::options_error`); a
/// reader exception is `BLAST query error: ...` (`ReadError::into_app_error`). The source's
/// `ReaderConfig::subject` carries `read_proteins` and the data loaders of `dlconfig`;
/// `-parse_deflines` is rejected by the programs and `gaps_to_Ns` is the mapper's; LOSAT keeps
/// the lowercase mask in every record (`FastaRecord::sequence`).
pub fn read_subjects<R: Read>(
    source: &mut FastaInputSource<R>,
    range: Option<&crate::blastinput::seq_range::SequenceRange>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
) -> anyhow::Result<Vec<FastaRecord>> {
    use crate::blastinput::seq_range::{record_interval, RecordInterval};
    let mut records = Vec::new();
    while !source.end() {
        match source.next_sequence(warn) {
            Ok(record) => {
                if let Some(range) = range {
                    if record_interval(range, record.sequence.len()) == RecordInterval::PastEnd {
                        return Err(crate::blastinput::app::options_error(
                            "Invalid from coordinate (greater than sequence length)",
                        ));
                    }
                }
                records.push(record);
            }
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => break,
            Err(error) => return Err(error.into_app_error()),
        }
    }
    Ok(records)
}

/// The queries of an input read to its end, as NCBI's batches read them one record at a
/// time (`CBlastInput::GetNextSeqBatch`): the records, and how the reading ended. The
/// caller cuts the batches (`query_batch.rs`) and writes each record's messages when its
/// batch is read.
#[derive(Debug)]
pub struct QueryRecords {
    pub records: Vec<FastaRecord>,
    pub end: QueryEnd,
}

/// How the reading of the queries ended.
#[derive(Debug)]
pub enum QueryEnd {
    /// `End()`: the line reader is at the end of the input.
    Input,
    /// `eEOF` from the reader (only blank and comment lines were left; no record follows):
    /// the batch being read ends there.
    EofError,
    /// Another reader error, thrown while reading the record after `records`, with the
    /// messages written before it.
    Error { error: ReadError, warnings: Vec<u8> },
}

/// Reads all queries (`GetNextSeqBatch` until `End()`): the messages are kept with the
/// records (written later, batch by batch, by the caller).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:134-171
/// ```c++
/// CRef<CBlastQueryVector>
/// CBlastInput::GetNextSeqBatch(CScope& scope)
/// {
///     CRef<CBlastQueryVector> retval(new CBlastQueryVector);
///     TSeqPos size_read = 0;
///
///     while (size_read < GetBatchSize()) {
///
///         if (End())
///             break;
///
///         CRef<CBlastSearchQuery> q;
///         try { q.Reset(m_Source->GetNextSequence(scope)); }
///         catch (const CObjReaderParseException& e) {
///             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
///                 break;
///             }
///             throw;
///         }
///         catch (const exception&) {
///             continue; //SB-2307. ignore well formed, not found accession
///         }
/// ```
/// The `catch (const exception&)` skips a query whose Seq-id was not found or whose
/// range does not fit (`seq_range.rs`); LOSAT rejects Seq-id lines.
pub fn read_queries<R: Read>(source: &mut FastaInputSource<R>) -> QueryRecords {
    let mut records = Vec::new();
    loop {
        if source.end() {
            return QueryRecords {
                records,
                end: QueryEnd::Input,
            };
        }
        let mut warnings = Vec::new();
        match source.next_sequence(&mut |message: &[u8]| {
            warnings.extend_from_slice(message);
            Ok(())
        }) {
            Ok(record) => records.push(record),
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => {
                return QueryRecords {
                    records,
                    end: QueryEnd::EofError,
                }
            }
            Err(error) => {
                return QueryRecords {
                    records,
                    end: QueryEnd::Error { error, warnings },
                }
            }
        }
    }
}
