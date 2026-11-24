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
use crate::seqid::accguide::SeqIdWarning;
use crate::seqid::ParsedSeqId;
use crate::seqmasks_io::bdb_reader::{prints_title, BdbError, MaskBdbReader};
use crate::seqmasks_io::fasta_reader::{MaskFastaReader, MaskReadError};
use crate::seqmasks_io::run::{
    finish, post_read_texts, serial_write_failure, write_failure, RunFailure,
};
use crate::seqmasks_io::writer::{MaskWriter, MaskedBioseq, Serial, WriteError};
use crate::sequence::title::get_title_without_source;
use crate::utils::winmask::{MaskerParams, PctSettings, SeqMasker, UnitCounts, WinMaskError};
use crate::writedb::fasta_source::reader::ReaderWarning;

use super::ids::IdLists;
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
    /// `sequence::GetTitle` of the Bioseq, which a list matched as text reads; made when the
    /// source was asked for it (`SequenceSource::make_titles`).
    pub get_title: Option<Vec<u8>>,
    pub iupac: Vec<u8>,
    pub texts: Vec<String>,
}

/// How a run stops early.
pub(super) enum Failure {
    /// An NCBI exception that leaves `Run`: reported by `NCBI_REPORT_EXCEPTION_X` (subcode 16,
    /// title "Application's execution failed"), exit status 3.
    Ncbi(NcbiException),
    /// A `std::exception` that leaves `Run`: posted by `ERR_POST_X` (subcode 7), exit status 3.
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
/// Writes the report of `failure` to the diagnostics; the exit status.
pub(super) fn report(failure: Failure, diag: &mut Diag) -> i32 {
    match failure {
        Failure::Ncbi(exception) => {
            diag.write_text(&exception.report_default_ex(
                106,
                16,
                "Application's execution failed",
            ));
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

/// `CMaskReader::GetNextSequence` of a source of sequences; `out` is the output, which comes
/// out before a post of the reader (NCBI's `cerr` is tied to `cout`).
pub(super) trait SequenceSource {
    fn next_record(
        &mut self,
        out: &mut Ostream<'_>,
        diag: &mut Diag,
    ) -> Result<Option<Record>, Failure>;
    /// From now on the records carry their titles (`Record::get_title`).
    fn make_titles(&mut self);
}

/// The reader of `-infmt` (`CWinMaskConfig`'s `reader`): `CMaskFastaReader` of the opened input
/// (nucleotide, with `-parse_seqids`), or `CMaskBDBReader` of the database `-in`; none for
/// `seqids` (`input` is `None` only for `blastdb`).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:384-397
pub(super) fn reader<'a, R: Read + 'a>(
    options: &Options<'_>,
    input: Option<R>,
) -> Result<Option<Box<dyn SequenceSource + 'a>>, Failure> {
    match (options.infmt, input) {
        ("fasta", Some(input)) => fasta_source(input, options.parse_seqids).map(Some),
        ("blastdb", _) => database_source(
            options.input,
            options.blastdb_registry,
            prints_title(options.outfmt),
        )
        .map(Some),
        _ => Ok(None),
    }
}

/// `CMaskFastaReader` of `input`, nucleotide, with `parse_seqids`: a failed first read is the stream's
/// `std::ios_base::failure`.
pub(super) fn fasta_source<'a, R: Read + 'a>(
    input: R,
    parse_seqids: bool,
) -> Result<Box<dyn SequenceSource + 'a>, Failure> {
    let reader = MaskFastaReader::new(input, false, parse_seqids).map_err(read_failure)?;
    Ok(Box::new(FastaSource {
        reader,
        titles: false,
    }))
}

/// `CMaskBDBReader` of the nucleotide database `name` on the search path of `registry`
/// (`Options::blastdb_registry`); `title_printed` as `MaskBdbReader::open`'s.
pub(super) fn database_source(
    name: &str,
    registry: Result<Option<&str>, &str>,
    title_printed: bool,
) -> Result<Box<dyn SequenceSource>, Failure> {
    Ok(Box::new(DatabaseSource::open(
        name,
        registry,
        title_printed,
    )?))
}

/// A failed read of `-in`: the `std::ios_base::failure` of NCBI's stream.
fn read_failure(error: std::io::Error) -> Failure {
    crate::seqmasks_io::run::read_failure(error, super::PROGRAM)
}

struct FastaSource<R: Read> {
    reader: MaskFastaReader<R>,
    /// `make_titles` was called.
    titles: bool,
}

impl<R: Read> SequenceSource for FastaSource<R> {
    // windowmasker posts at the level Warning (`SetDiagPostLevel` to `eDiag_Warning`,
    // win_mask_app.cpp:118): the reader's warnings are written, each on its line (the line
    // warnings, `LOG_POST`) or as `ERR_POST_X` of the Seq-id parser (error code 1306). The
    // output written so far comes out first, whatever its last byte (`cerr` is tied to `cout`);
    // a closed pipe ends the run there, as SIGPIPE ends NCBI's.
    fn next_record(
        &mut self,
        out: &mut Ostream<'_>,
        diag: &mut Diag,
    ) -> Result<Option<Record>, Failure> {
        let mut closed = false;
        let record = self.reader.next(&mut |warning| {
            let _ = out.flush();
            closed |= out.closed();
            if closed {
                return;
            }
            match warning {
                ReaderWarning::Line(line) => diag.post_bytes(line),
                ReaderWarning::SeqId(warning) => post_seq_id_warning(diag, &warning),
            }
        });
        if closed {
            return Err(Failure::closed_pipe());
        }
        let record = record.map_err(|error| match error {
            MaskReadError::Ncbi(exception) => Failure::Ncbi(exception),
            MaskReadError::Read(error) => read_failure(error),
            MaskReadError::OutOfMemory => Failure::Std("std::bad_alloc".to_string()),
            MaskReadError::Std(what) => Failure::Std(what),
            MaskReadError::Refused(what) => Failure::Refused(unsupported(&what)),
            MaskReadError::SeqMap { error, .. } => {
                Failure::Ncbi(crate::seqmasks_io::fasta_reader::seq_map_exception(error))
            }
        })?;
        let titles = self.titles;
        Ok(record.map(|record| Record {
            get_title: titles.then(|| {
                get_title_without_source(&record.ids, record.title.as_deref().unwrap_or_default())
            }),
            ids: record.ids,
            title: record.title,
            taxname: None,
            has_source: false,
            iupac: record.letters,
            texts: Vec::new(),
        }))
    }

    fn make_titles(&mut self) {
        self.titles = true;
    }
}

/// A warning of NCBI's Seq-id code (`ERR_POST_X` of error code 1306, at windowmasker's post
/// level Warning) as windowmasker posts it: `Warning: (1306.<subcode>) <message>`, no prefix.
pub(super) fn post_seq_id_warning(diag: &mut Diag, warning: &SeqIdWarning) {
    let mut line = format!("Warning: (1306.{}) ", warning.subcode).into_bytes();
    line.extend_from_slice(&warning.message);
    diag.post_bytes(&line)
}

/// `CMaskBDBReader` of `name`, a nucleotide database (`seqmasks_io::bdb_reader`).
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
    fn next_record(
        &mut self,
        _out: &mut Ostream<'_>,
        _diag: &mut Diag,
    ) -> Result<Option<Record>, Failure> {
        let record = self.reader.next_record().map_err(database_failure)?;
        Ok(record.map(|record| Record {
            ids: record.ids,
            title: Some(record.title),
            taxname: record.taxname,
            has_source: record.has_source,
            get_title: record.get_title,
            iupac: record.residues,
            texts: record.texts,
        }))
    }

    fn make_titles(&mut self) {
        self.reader.ask_get_title();
    }
}

/// `CWinMaskConfig::x_GetWriter`: the writer of `-outfmt`.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:282-318
/// (`-outfmt` takes only the eight values of its constraint, of which none starts with
/// `seqloc_asn1_binary` or `maskinfo_asn1_binary`. `BuildAlgorithmParametersString` gives an
/// empty string: windowmasker has none of `locut`, `hicut`, `window` with `level` and
/// `linker`, mask_writer_blastdb_maskinfo.cpp:180-202.)
pub(super) fn writer(outfmt: &str) -> MaskWriter {
    let mask_info =
        |serial| MaskWriter::mask_info(Some(serial), 3, FilterProgram::Windowmasker, Vec::new());
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
fn masker(options: &Options<'_>, diag: &mut Diag) -> Result<SeqMasker, Failure> {
    let file = format!("-ustat {}", options.ustat);
    let counts = open_counts(options.ustat, &file, &options.params, &options.pct, diag)?;
    SeqMasker::new(Arc::new(counts), &options.params).map_err(|error| winmask_failure(&file, error))
}

/// `CSeqMaskerIstatFactory::create` of the counts file `name` with the thresholds `params` and the
/// percentages `pct`, its warnings posted as NCBI posts them: the factory's at once, the
/// readers' after. `file` names the file in LOSAT's refusals (`-ustat <name>`).
pub(super) fn open_counts(
    name: &str,
    file: &str,
    params: &MaskerParams,
    pct: &PctSettings,
    diag: &mut Diag,
) -> Result<UnitCounts, Failure> {
    let counts = UnitCounts::open_with_pct(Path::new(name), params, pct, &mut |warning| {
        diag.warning_post(&warning)
    });
    let counts = match counts {
        Ok(counts) => counts,
        Err(error) => {
            for warning in error.warnings() {
                diag.warning_post(warning);
            }
            return Err(winmask_failure(file, error));
        }
    };
    for warning in counts.warnings() {
        diag.warning_post(warning);
    }
    Ok(counts)
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_hash.cpp:94-100
/// (`CSeqMaskerUsetHash::Exception::GetErrCodeString`: "bad index", :43-50.)
const BAD_INDEX: ThrowSite = ThrowSite {
    file: "src/algo/winmask/seq_masker_uset_hash.cpp",
    line: 99,
    kind: "Exception::bad index",
    function: Some("ncbi::CSeqMaskerUsetHash::get_info"),
};

/// An error of the counts file `file` (`-ustat <name>`, `-in <name>`) or of the masker: NCBI's
/// exception, or LOSAT's refusal where NCBI's result is undefined.
pub(super) fn winmask_failure(file: &str, error: WinMaskError) -> Failure {
    if let Some(exception) = error.create_exception() {
        return Failure::Ncbi(exception);
    }
    match error {
        WinMaskError::BadIndex { .. } => {
            Failure::Ncbi(NcbiException::new(BAD_INDEX, error.to_string()))
        }
        WinMaskError::PctLine { .. } => Failure::Refused(unsupported(&format!("{file}: {error};"))),
        _ => Failure::Refused(unsupported(&format!(
            "{file}, which gives an error that NCBI BLAST+ does not report as such ({error}),"
        ))),
    }
}

// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:238-280
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
    while let Some(record) = source.next_record(out, diag)? {
        post_read_texts::<Failure>(&record.texts, out, diag)?;
        if record.iupac.is_empty() || !lists.consider(&record) {
            continue;
        }
        let mut mask = masker
            .mask(&record.iupac)
            .map_err(|error| winmask_failure(&format!("-ustat {}", options.ustat), error))?;
        if let Some(duster) = duster.as_mut() {
            let dust = duster.mask(&record.iupac, &mask);
            merge_mask_info(&mut mask, &dust);
        }
        let seq = MaskedBioseq {
            ids: &record.ids,
            handle_seq_id: lists.handle_seq_id(&record.ids),
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
