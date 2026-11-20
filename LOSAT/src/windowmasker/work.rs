//! `CWinMaskApplication::Run` with a counts file (`-ustat`): the reader of `-infmt`, the writer
//! of `-outfmt`, the id lists, and the masks of each nucleotide sequence (WindowMasker, then
//! DUST with `-dust T`).

use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;

use crate::ncbi_app::files::Ostream;
use crate::ncbi_app::Diag;
use crate::ncbi_exception::{err_post_x, NcbiException, ThrowSite};
use crate::seqdb::mask::FilterProgram;
use crate::seqdb::SeqType;
use crate::seqid::ParsedSeqId;
use crate::seqmasks_io::bdb_reader::{prints_title, BdbError, MaskBdbReader};
use crate::seqmasks_io::fasta_reader::{MaskFastaReader, MaskReadError};
use crate::seqmasks_io::run::{
    finish, post_read_texts, serial_write_failure, write_failure, RunFailure,
};
use crate::seqmasks_io::writer::{MaskWriter, MaskedBioseq, Serial, WriteError};
use crate::utils::winmask::{MaskerParams, PctSettings, SeqMasker, UnitCounts, WinMaskError};
use crate::writedb::fasta_source::reader::ReaderWarning;

use super::ids::{get_title, IdLists, Title};
use super::sdust::{merge_mask_info, SDustMasker};

/// The settings of the masking run, as `CWinMaskConfig` reads them from NCBI's `CArgs`.
pub(super) struct Options<'a> {
    pub ustat: &'a str,
    pub input: &'a str,
    pub infmt: &'a str,
    pub outfmt: &'a str,
    pub parse_seqids: bool,
    pub params: MaskerParams,
    pub pct: PctSettings,
    /// `-dust_level` when `-dust` is true (`eGenerateMasksWithDuster`).
    pub dust_level: Option<i32>,
    /// The `BLASTDB` of the configuration (`ApplicationSettings::blastdb`).
    pub blastdb_registry: Result<Option<&'a str>, &'a str>,
}

/// One sequence of the input (`CBioseq`).
pub(super) struct Record {
    pub ids: Vec<ParsedSeqId>,
    pub title: Option<Vec<u8>>,
    pub taxname: Option<Arc<[u8]>>,
    pub has_source: bool,
    /// Whether the Bioseq has a source (`BdbRecord::any_source`), for `GetTitle`.
    pub any_source: bool,
    pub iupac: Vec<u8>,
    pub texts: Vec<String>,
}

/// How a run stops early.
pub(super) enum Failure {
    /// An NCBI exception that leaves `Run`: `NCBI_REPORT_EXCEPTION_X(16, "Application's
    /// execution failed", e)`, exit status 3.
    Ncbi(NcbiException),
    /// A `std::exception` that leaves `Run`: `ERR_POST_X(7, ...)`, exit status 3.
    Std(String),
    /// A case LOSAT does not reproduce.
    Refused(String),
    /// A `CArgException` thrown in `Run` (`AsOutputFile`).
    Arg(crate::ncbi_app::args::ArgException),
    /// A write found the pipe closed where SIGPIPE ends NCBI's process.
    ClosedPipe,
}

impl From<WriteError> for Failure {
    fn from(error: WriteError) -> Self {
        write_failure(error, super::PROGRAM)
    }
}

impl RunFailure for Failure {
    fn refused(message: String) -> Self {
        Failure::Refused(message)
    }

    fn closed_pipe() -> Self {
        Failure::ClosedPipe
    }

    fn std_exception(what: String) -> Self {
        Failure::Std(what)
    }

    fn stops_writing(&self) -> bool {
        matches!(self, Failure::Refused(_) | Failure::ClosedPipe)
    }
}

/// LOSAT's refusal of `what` (a noun phrase).
pub(super) fn unsupported(what: &str) -> String {
    crate::seqmasks_io::run::unsupported(what, super::PROGRAM)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:668-690 (AppMain's handlers of `Run`)
// ```c++
//         catch (CException& e) {
//             NCBI_REPORT_EXCEPTION_X(16,
//                                     "Application's execution failed", e);
//             *got_exception = true;
//             *exit_code = 3;
//         }
//         catch (exception& e) {
//             ERR_POST_X(7, "Application's execution failed: " << e.what());
//             *got_exception = true;
//             *exit_code = 3;
//         }
// ```
/// Writes the report of `failure` to the diagnostics; the exit status.
pub(super) fn report(failure: Failure, diag: &mut Diag) -> i32 {
    match failure {
        Failure::Ncbi(exception) => {
            diag.write_text(
                exception
                    .report_default_ex(106, 16, "Application's execution failed")
                    .as_bytes(),
            );
            3
        }
        Failure::Std(what) => {
            diag.write_text(
                err_post_x(106, 7, &format!("Application's execution failed: {what}")).as_bytes(),
            );
            3
        }
        Failure::Refused(reason) => diag.refuse(&reason),
        Failure::Arg(error) => diag.run_arg_exception(error),
        Failure::ClosedPipe => 0,
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:486-506
/// ```c++
/// const char * CWinMaskConfigException::GetErrCodeString() const
/// {
///     switch( GetErrCode() )
///     {
///     case eInputOpenFail:
///
///         return "can not open input stream";
/// ...
///     case eInconsistentOptions:
///
///         return "inconsistent program options";
/// ```
/// `CWinMaskConfigException` thrown at `line` of win_mask_config.cpp in `function`.
pub(super) fn config_exception(
    line: u32,
    kind: &'static str,
    function: &'static str,
    message: &str,
) -> NcbiException {
    NcbiException::new(
        ThrowSite {
            file: "src/algo/winmask/win_mask_config.cpp",
            line,
            kind,
            function: Some(function),
        },
        message,
    )
}

/// The `kind` of `config_exception` for `eInconsistentOptions`.
pub(super) const INCONSISTENT: &str = "CWinMaskConfigException::inconsistent program options";

/// The function of the `CWinMaskConfig` constructor in NCBI's reports.
pub(super) const CONFIG: &str = "ncbi::CWinMaskConfig::CWinMaskConfig";

/// `CMaskReader::GetNextSequence` of a source of sequences.
pub(super) trait SequenceSource {
    fn next_record(&mut self, diag: &mut Diag) -> Result<Option<Record>, Failure>;
}

/// The reader of `-infmt` (`CWinMaskConfig`'s `reader`): `CMaskFastaReader( *is, true,
/// args["parse_seqids"] )` of the opened input, or `CMaskBDBReader( args[kInput].AsString() )`;
/// none for `seqids` (`input` is `None` only for `blastdb`).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:384-397
/// ```c++
///         if(determine_input && iformatstr != "seqids"){
///             if( iformatstr == "fasta" )
///                 reader = new CMaskFastaReader( *is, true, args["parse_seqids"] );
///             else if( iformatstr == "blastdb" )
///                 reader = new CMaskBDBReader( args[kInput].AsString() );
/// ```
pub(super) fn reader<'a, R: Read + 'a>(
    options: &Options<'_>,
    input: Option<R>,
) -> Result<Option<Box<dyn SequenceSource + 'a>>, Failure> {
    match (options.infmt, input) {
        ("fasta", Some(input)) => {
            let reader = MaskFastaReader::new(input, false, options.parse_seqids)
                .map_err(read_failure)?;
            Ok(Some(Box::new(FastaSource { reader })))
        }
        ("blastdb", _) => Ok(Some(Box::new(DatabaseSource::open(
            options.input,
            options.blastdb_registry,
            prints_title(options.outfmt),
        )?))),
        _ => Ok(None),
    }
}

/// A failed read of `-in`: the `std::ios_base::failure` of NCBI's stream.
fn read_failure(error: std::io::Error) -> Failure {
    crate::seqmasks_io::run::read_failure(error, super::PROGRAM)
}

struct FastaSource<R: Read> {
    reader: MaskFastaReader<R>,
}

impl<R: Read> SequenceSource for FastaSource<R> {
    // windowmasker posts at the level Warning (`SetDiagPostLevel( eDiag_Warning )`,
    // win_mask_app.cpp:118): the reader's warnings are written, each on its line (the line
    // warnings, `LOG_POST`) or as `ERR_POST_X` of the Seq-id parser (error code 1306).
    fn next_record(&mut self, diag: &mut Diag) -> Result<Option<Record>, Failure> {
        let record = self
            .reader
            .next(&mut |warning| match warning {
                ReaderWarning::Line(line) => diag.post_bytes(line),
                ReaderWarning::SeqId(warning) => diag.post(&format!(
                    "Warning: (1306.{}) {}",
                    warning.subcode, warning.message
                )),
            })
            .map_err(|error| match error {
                MaskReadError::Ncbi(exception) => Failure::Ncbi(exception),
                MaskReadError::Read(error) => read_failure(error),
                MaskReadError::OutOfMemory => Failure::Std("std::bad_alloc".to_string()),
                MaskReadError::Refused(what) => Failure::Refused(unsupported(&what)),
            })?;
        Ok(record.map(|record| Record {
            ids: record.ids,
            title: record.title,
            taxname: None,
            has_source: false,
            any_source: false,
            iupac: record.letters,
            texts: Vec::new(),
        }))
    }
}

/// `CMaskBDBReader(name)` of a nucleotide database (`seqmasks_io::bdb_reader`).
struct DatabaseSource {
    reader: MaskBdbReader,
}

impl DatabaseSource {
    fn open(
        name: &str,
        registry: Result<Option<&str>, &str>,
        title_printed: bool,
    ) -> Result<Self, Failure> {
        let reader = MaskBdbReader::open(
            name,
            SeqType::Nucleotide,
            registry,
            title_printed,
            super::PROGRAM,
        )
        .map_err(database_failure)?;
        Ok(DatabaseSource { reader })
    }
}

/// An error of the database reader, reported as `AppMain` reports an exception of `Run`, or
/// LOSAT's refusal.
fn database_failure(error: BdbError) -> Failure {
    match error {
        BdbError::Ncbi(exception) => Failure::Ncbi(exception),
        BdbError::Refused(reason) => Failure::Refused(reason),
        BdbError::Unsupported(what) => Failure::Refused(unsupported(&what)),
    }
}

impl SequenceSource for DatabaseSource {
    fn next_record(&mut self, _diag: &mut Diag) -> Result<Option<Record>, Failure> {
        let record = self.reader.next_record().map_err(database_failure)?;
        Ok(record.map(|record| Record {
            ids: record.ids,
            title: Some(record.title),
            taxname: record.taxname,
            has_source: record.has_source,
            any_source: record.any_source,
            iupac: record.residues,
            texts: record.texts,
        }))
    }
}

/// `CWinMaskConfig::x_GetWriter`: the writer of `-outfmt`.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:282-318
/// ```c++
///     if (format == "interval") {
///         CNcbiOstream& output = args[kOutput].AsOutputFile();
///         retval = new CMaskWriterInt(output);
///     } else if (format == "fasta") {
///         CNcbiOstream& output = args[kOutput].AsOutputFile();
///         retval = new CMaskWriterFasta(output);
///     } else if (NStr::StartsWith(format, "seqloc_asn1_binary")) {
///         ...
///     } else if (NStr::StartsWith(format, "seqloc_")) {
///         CNcbiOstream& output = args[kOutput].AsOutputFile();
///         retval = new CMaskWriterSeqLoc(output, format);
///     } else if (NStr::StartsWith(format, "maskinfo_asn1_binary")) {
///         ...
///     } else if (NStr::StartsWith(format, "maskinfo_")) {
///         CNcbiOstream& output = args[kOutput].AsOutputFile();
///         retval =
///             new CMaskWriterBlastDbMaskInfo(output, format, 3,
///                                eBlast_filter_program_windowmasker,
///                                BuildAlgorithmParametersString(args));
/// ```
/// (`-outfmt` takes only the eight values of its constraint, of which none starts with
/// `seqloc_asn1_binary` or `maskinfo_asn1_binary`. `BuildAlgorithmParametersString` gives an
/// empty string: windowmasker has none of `locut`, `hicut`, `window` with `level` and
/// `linker`, mask_writer_blastdb_maskinfo.cpp:180-202.)
pub(super) fn writer(outfmt: &str) -> MaskWriter {
    let mask_info = |serial| {
        MaskWriter::mask_info(Some(serial), 3, FilterProgram::Windowmasker, Vec::new())
    };
    match outfmt {
        "interval" => MaskWriter::Interval,
        "fasta" => MaskWriter::Fasta {
            hard_masking: false,
        },
        "seqloc_asn1_bin" => MaskWriter::SeqLoc(Serial::AsnBinary),
        "seqloc_asn1_text" => MaskWriter::SeqLoc(Serial::AsnText),
        "seqloc_xml" => MaskWriter::SeqLoc(Serial::Xml),
        "maskinfo_asn1_bin" => mask_info(Serial::AsnBinary),
        "maskinfo_asn1_text" => mask_info(Serial::AsnText),
        "maskinfo_xml" => mask_info(Serial::Xml),
        other => unreachable!("-outfmt {other} passed its constraint"),
    }
}

/// `CWinMaskApplication::Run` after the configuration: `source` is the reader (none for
/// `seqids`), `writer` the writer, `lists` the id lists; the exit status. The `maskinfo_*`
/// writer writes its object also when an exception leaves `Run` (the destructor of the
/// configuration deletes the writer).
///
/// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:203-206
/// ```c++
///     if(aConfig.InFmt() == "seqids"){
///         ERR_POST(Error << "windowmasker with seqids input not implemented yet");
///         return 1;
///     }
/// ```
pub(super) fn run(
    options: &Options<'_>,
    source: Option<Box<dyn SequenceSource + '_>>,
    writer: MaskWriter,
    lists: &IdLists,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> Result<i32, Failure> {
    let mut out = Ostream::new(out);
    let mut writer = writer;
    let (result, status) = match source {
        None => {
            diag.err_post("windowmasker with seqids input not implemented yet");
            (Ok(()), 1)
        }
        Some(mut source) => (
            mask_input(options, &mut *source, &mut writer, lists, &mut out, diag),
            0,
        ),
    };
    finish(
        writer,
        &mut out,
        diag,
        result,
        options.outfmt,
        super::PROGRAM,
    )?;
    Ok(status)
}

/// The masker of the run: `CSeqMasker`'s counts (`CSeqMaskerIstatFactory::create`, with the
/// `-t_*_pct` settings) and its window check. The warnings of the counts file are posted as
/// NCBI posts them: the factory's at once, the readers' after.
///
/// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:208-232
/// ```c++
///     CMaskReader & theReader = aConfig.Reader();
///     CMaskWriter & theWriter = aConfig.Writer();
///     CSeqMasker theMasker( aConfig.LStatName(),
///                           aConfig.WindowSize(),
///                           ...
///                           aConfig.MinScorePct(),
///                           aConfig.ExtendScorePct(),
///                           aConfig.ThresScorePct(),
///                           aConfig.MaxScorePct() );
/// ```
fn masker(options: &Options<'_>, diag: &mut Diag) -> Result<SeqMasker, Failure> {
    let counts = UnitCounts::open_with_pct(
        Path::new(options.ustat),
        &options.params,
        &options.pct,
        &mut |warning| diag.warning_post(&warning),
    );
    let counts = match counts {
        Ok(counts) => counts,
        Err(error) => {
            for warning in error.warnings() {
                diag.warning_post(warning);
            }
            return Err(winmask_failure(options.ustat, error));
        }
    };
    for warning in counts.warnings() {
        diag.warning_post(warning);
    }
    SeqMasker::new(Arc::new(counts), &options.params)
        .map_err(|error| winmask_failure(options.ustat, error))
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_hash.cpp:94-100
/// ```c++
///         if( (hval>>bc) + coll > M )
///         {
///             ostringstream r;
///             r << "bad index at key " << hash.first
///               << " : " << htp[hash.first];
///             NCBI_THROW( Exception, eBadIndex, r.str() );
///         }
/// ```
/// (`CSeqMaskerUsetHash::Exception::GetErrCodeString`: "bad index", :43-50.)
const BAD_INDEX: ThrowSite = ThrowSite {
    file: "src/algo/winmask/seq_masker_uset_hash.cpp",
    line: 99,
    kind: "Exception::bad index",
    function: Some("ncbi::CSeqMaskerUsetHash::get_info"),
};

/// An error of the counts file `ustat` or of the masker: NCBI's exception, or LOSAT's refusal
/// where NCBI's result is undefined.
fn winmask_failure(ustat: &str, error: WinMaskError) -> Failure {
    if let Some(exception) = error.create_exception() {
        return Failure::Ncbi(exception);
    }
    match error {
        WinMaskError::BadIndex { .. } => {
            Failure::Ncbi(NcbiException::new(BAD_INDEX, error.to_string()))
        }
        WinMaskError::PctLine { .. } => Failure::Refused(unsupported(&format!(
            "-ustat {ustat} with a -t_*_pct setting: {error};"
        ))),
        _ => Failure::Refused(unsupported(&format!(
            "-ustat {ustat}, which gives an error that NCBI BLAST+ does not report as such ({error}),"
        ))),
    }
}

// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:238-280
// ```c++
//     if( aConfig.AppType() == CWinMaskConfig::eGenerateMasksWithDuster )
//         duster = new CSDustMasker( aConfig.DustWindow(),
//                                    aConfig.DustLevel(),
//                                    aConfig.DustLinker() );
//
//     while( (aSeqEntry = theReader.GetNextSequence()).NotEmpty() )
//     {
//         if( aSeqEntry->Which() == CSeq_entry::e_not_set ) continue;
//         ...
//         for ( ;  bs_iter;  ++bs_iter) {
//             CBioseq_Handle bsh = *bs_iter;
//             if (bsh.GetBioseqLength() == 0) {
//                 continue;
//             }
//
//             if( CWinMaskUtil::consider( bsh, ids, exclude_ids ) )
//             {
//                 ...
//                 unique_ptr< CSeqMasker::TMaskList > mask_info( theMasker( data ) );
//                 CSeqMasker::TMaskList dummy;
//
//                 if( duster != 0 ) // Dust and merge with mask_info
//                 {
//                     unique_ptr< CSeqMasker::TMaskList > dust_info(
//                         (*duster)( data, *mask_info.get() ) );
//                     CSeqMasker::MergeMaskInfo( mask_info.get(), dust_info.get() );
//                 }
//
//                 // theWriter.Print( bsh, *mask_info, aConfig.MatchId() );
//                 theWriter.Print( bsh, *mask_info, GetArgs()["parse_seqids"] );
// ```
/// The masks of every sequence of `source`, printed by `writer`.
fn mask_input(
    options: &Options<'_>,
    source: &mut (dyn SequenceSource + '_),
    writer: &mut MaskWriter,
    lists: &IdLists,
    out: &mut Ostream<'_>,
    diag: &mut Diag,
) -> Result<(), Failure> {
    let masker = masker(options, diag)?;
    let mut duster = options.dust_level.map(SDustMasker::new);
    while let Some(record) = source.next_record(diag)? {
        post_read_texts::<Failure>(&record.texts, out, diag)?;
        if record.iupac.is_empty() || !considered(lists, &record)? {
            continue;
        }
        let mut mask = masker
            .mask(&record.iupac)
            .map_err(|error| winmask_failure(options.ustat, error))?;
        if let Some(duster) = duster.as_mut() {
            let dust = duster.mask(&record.iupac, &mask);
            merge_mask_info(&mut mask, &dust);
        }
        let seq = MaskedBioseq {
            ids: &record.ids,
            title: record.title.as_deref(),
            taxname: record.taxname.as_ref(),
            has_source: record.has_source,
            iupac: &record.iupac,
            protein: false,
        };
        writer.print(&seq, &mask, options.parse_seqids, out, diag)?;
        serial_write_failure::<Failure>(out, options.outfmt, super::PROGRAM)?;
    }
    Ok(())
}

/// `CWinMaskUtil::consider` of `record`, which asks for its title only when a list is used.
fn considered(lists: &IdLists, record: &Record) -> Result<bool, Failure> {
    if !lists.active() {
        return Ok(true);
    }
    match get_title(&record.ids, record.title.as_deref(), record.any_source) {
        Title::Text(title) => Ok(lists.consider(&title)),
        Title::FromSource => Err(Failure::Refused(unsupported(
            "-ids or -exclude_ids with a database sequence without a title whose definition line has a tax id (NCBI matches the ids against a title made from the taxonomy names)",
        ))),
    }
}
