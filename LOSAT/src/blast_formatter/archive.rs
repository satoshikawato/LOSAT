//! blast_formatter's run over the objects of a BLAST archive (`-archive`): NCBI's
//! `CBlastFormatterApp::Run` with `CRemoteBlast::LoadFromArchive` and `PrintFormattedOutput`
//! (c++/src/app/blast/blast_formatter.cpp:329-531) over the archive reader
//! (`report::archive::read`): the formats 8, 9 and 12, which write each query's Seq-aligns as
//! the archive holds them, cut to the hit list size, inside the one Seq-annot wrapper of the
//! search programs (`report::serial::seqalign::write_document_to`, decision F3-2 of line RG), and
//! the format 11, which writes the archive again with the search programs' writer
//! (`archive_out`, `report::archive::writer::write_archive_document`).
//!
//! The object loop keeps what NCBI's `CRemoteBlast` keeps from one object to the next: the
//! options of the first search object (`GetSearchOptions` builds them once: the program, the
//! taxids, the database's filtering algorithm), the hit list size that `ExtractAlgorithmOptions`
//! set last, the database that the last database object named (`m_Dbs`, which an object of
//! `-subject` sequences does not clear) and the `-out` stream, opened when the first object is
//! formatted. Where NCBI's result would depend on its servers (a database that cannot be opened
//! here), on memory its input does not decide, or on a form LOSAT does not read, the run stops
//! with LOSAT's refusal (`refusals::Refusal`) at the point where NCBI meets it.

use std::cell::RefCell;
use std::ffi::OsStr;
use std::io::{Read, Write};
use std::path::Path;

use super::archive_out;
use super::formatting::{extract_formatting, sam_check};
use super::refusals::{self, Refusal};
use super::search_program::{blast_type, program_type};
use super::{Options, PROGRAM};
use crate::algorithm::common::database::{
    database_error, database_error_of, list_warning_line, open_search_database, pairwise_database,
};
use crate::blastinput::app::{
    engine_error, options_error_bytes, output_format, prefetch_seqs_limit, xinclude_check_for,
};
use crate::blastinput::db_mask::name_mask_error;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::{command_line_argument_error, NativeError};
use crate::ncbi_app::files::{self, OutputDevice};
use crate::ncbi_app::Diag;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::archive::fields;
use crate::report::archive::options::{ParamValue, RemoteOptions};
use crate::report::archive::read::options::{archive_request, BuiltOptions};
use crate::report::archive::read::request::{ParamSet, Queries, Subject};
use crate::report::archive::read::{
    seq_align_id, Archive, ArchiveReader, Message, Search, Severity, Stop,
};
use crate::report::archive::writer::FLUSH_BYTES;
use crate::report::run_files::{OpenFiles, OpenedFile};
use crate::report::serial::asn_input::Names;
use crate::report::serial::asn_types::read_xml::XmlStyle;
use crate::report::serial::asn_types::stream::{guess_encoding, Encoding, GuessStop};
use crate::report::serial::asn_types::table::member;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::write::write_value;
use crate::report::serial::asn_types::ReadError;
use crate::report::serial::encoder::{Encoder, SerialError};
use crate::report::serial::seqalign::{
    write_document, write_document_posting, AnnotDesc, DatabaseDesc, SeqAlignError,
    SeqAnnotDocument, SerialFormat,
};
use crate::seqdb::{DbList, SeqDb, SeqDbError, SeqType, TaxidList};
use crate::seqid::defline_set::{compare, Sic};
use crate::seqid::ParsedSeqId;
use crate::writedb::cout::CoutBuf;
use crate::writedb::metadata::{OStreamBuffer, ObjectOStream};
use crate::writedb::ofstream::{BasicStream, FileBuf, IoState, StreamBuf};

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:506-529
/// The formatting of the archive `archive` with the command line's `options` and the
/// application's `settings` (the registry's `[BLAST] BLASTDB`), its reports written to `-out`
/// (`out` is the standard output, `files` the run's files), its diagnostics on `diag`; the exit
/// status. The errors are
/// `CATCH_ALL`'s (`super::catch_all`), after what was written before them; the stream writes what
/// it holds when the run ends (`CArgs` closes the file of `-out`; `cout` is flushed at the
/// program's end).
pub(crate) fn format_archive(
    options: &Options,
    settings: &ApplicationSettings,
    archive: &mut dyn Read,
    out: &mut dyn Write,
    files: &OpenFiles,
    diag: &mut Diag,
) -> i32 {
    let mut run = Run {
        options,
        settings,
        output: Output {
            stdout: Some(out),
            files,
            stream: None,
            formatting: false,
        },
        search_options: None,
        algorithm_options: None,
        hitlist: 0,
        database: None,
        psi_step: None,
        opened: None,
        buf: Vec::new(),
    };
    match run.objects(archive, diag) {
        Ok(()) => 0,
        Err(error) => super::catch_all(diag, error),
    }
}

/// What the loop keeps from one archive object to the next (`CRemoteBlast`'s members, and the
/// `-out` stream of `CArgs`).
struct Run<'r, 'o> {
    options: &'r Options,
    settings: &'r ApplicationSettings,
    output: Output<'o>,
    /// The options of the first search object (`m_CBOH`, `m_TaxidList`,
    /// `m_DbFilteringAlgorithmId`, `m_Task`), which every later object is formatted with.
    search_options: Option<BuiltOptions>,
    /// The remote options of `m_CBOH` (`-outfmt 11`), built once on the first search object
    /// (`RemoteOptions::built`, `GetSearchOptions`) and set across objects (`HitlistSize`, RPS
    /// `MatrixName`).
    algorithm_options: Option<RemoteOptions>,
    /// The options' hit list size: the first object's, then what `ExtractAlgorithmOptions` set
    /// for the last formatted object (`opts.SetHitlistSize` on the kept options).
    hitlist: i32,
    /// `m_Dbs`: the name and residue type of the database of the last database object.
    database: Option<(Vec<u8>, bool)>,
    /// `Web_StepNumber` of `m_FormatOpts` (`remote_blast.cpp:1733-1737`), updated whenever a
    /// search object has `service == "psi"`; `Err(())` when the value is not an `integer`.
    psi_step: Option<Result<i32, ()>>,
    /// The database opened for an earlier object (`opened_database`).
    opened: Option<OpenedDatabase>,
    /// The encoder's buffer of one output document.
    buf: Vec<u8>,
}

/// The `-out` stream: the standard output until `-out` is opened (`AsOutputFile`, at the first
/// formatted object), then `cout` or the file's `ofstream` as the release's C++ runtime has them
/// (`writedb::ofstream`, `writedb::cout`: what they hold, when a write fails, their state and its
/// exceptions).
pub(super) struct Output<'o> {
    /// The standard output until `-out` is opened.
    stdout: Option<&'o mut dyn Write>,
    /// The run's files, where `-out` is opened (`OpenFiles::create`).
    files: &'o OpenFiles,
    /// The opened stream.
    stream: Option<BasicStream<OutBuf<'o>>>,
    /// Whether an object's formatter holds the stream (`begin_object`, `end_object`).
    formatting: bool,
}

/// The buffer of the `-out` stream: the file's (`writedb::ofstream::FileBuf`, the `ofstream` of
/// an `eOutputFile` argument) or the standard output's (`writedb::cout::CoutBuf`, NCBI's `cout`:
/// the C library's `stdout`, whose block is the device's, 4096 bytes on a pipe or `/dev/full`;
/// observed with blast_formatter 2.17.0, a Seq-annot of 2390 bytes failing at its close and one
/// of 4845 bytes in its write, line RG item F3, agent C6, 2026-11-23,
/// `remaining-gaps-20261112/F/f3/tmp/c6/pc`). Its stream can be set to throw (`exceptions`).
enum OutBuf<'o> {
    File(FileBuf<OpenedFile>),
    Cout(CoutBuf<'o>),
}

impl StreamBuf for OutBuf<'_> {
    const MASKABLE: bool = true;

    fn put_bytes(&mut self, data: &[u8]) -> bool {
        match self {
            OutBuf::File(buf) => StreamBuf::put_bytes(buf, data),
            OutBuf::Cout(buf) => StreamBuf::put_bytes(buf, data),
        }
    }

    fn put_byte(&mut self, byte: u8) -> bool {
        match self {
            OutBuf::File(buf) => StreamBuf::put_byte(buf, byte),
            OutBuf::Cout(buf) => StreamBuf::put_byte(buf, byte),
        }
    }

    fn sync(&mut self) -> bool {
        match self {
            OutBuf::File(buf) => StreamBuf::sync(buf),
            OutBuf::Cout(buf) => StreamBuf::sync(buf),
        }
    }

    /// Never: blast_formatter links NCBI's connection library, which masks SIGPIPE
    /// (ncbi_socket.c:865-873), so a write to a closed pipe is a failed write like the others
    /// (observed: blast_formatter 2.17.0 aborts on a closed standard output as on `/dev/full`;
    /// line RG item F3, agent C6, 2026-11-23).
    fn closed_pipe(&self) -> bool {
        false
    }
}

impl Drop for OutBuf<'_> {
    /// At the program's end the C library writes what `stdout` holds, whatever `cout`'s state
    /// (the file of `-out` is closed by its buffer's own destructor).
    fn drop(&mut self) {
        if let OutBuf::Cout(buf) = self {
            StreamBuf::sync(buf);
        }
    }
}

impl<'o> Output<'o> {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:335
    /// `args[kArgOutput].AsOutputFile()` of `-out` `value`: opened the first time (NCBI's
    /// `Argument "out". File is not accessible` there), the same stream after.
    fn open(&mut self, value: &OsStr) -> anyhow::Result<()> {
        if self.stream.is_some() {
            return Ok(());
        }
        let stdout = self
            .stdout
            .take()
            .expect("the standard output until -out is opened");
        let file = files::open_output_device("out", value, stdout, self.files)
            .map_err(|error| command_line_argument_error(error.message()))?;
        let buf = match file {
            OutputDevice::Stdout(out) => {
                OutBuf::Cout(CoutBuf::new(out, crate::cli::NcbiStdio::cout()))
            }
            OutputDevice::File(file) => OutBuf::File(FileBuf::of_file(file)),
            // LOSAT Web's refusal of a name that it does not write, in blastdbcmd's words
            // (`OutputDevice::Refused`); the command line and the library open every name.
            OutputDevice::Refused(reason) => {
                return Err(anyhow::anyhow!(crate::seqmasks_io::run::unsupported(
                    &reason, PROGRAM
                )))
            }
        };
        self.stream = Some(BasicStream::with_buf(buf));
        Ok(())
    }

    /// The opened stream and whether it is `cout`.
    fn opened(&mut self) -> (&mut BasicStream<OutBuf<'o>>, bool) {
        let stream = self
            .stream
            .as_mut()
            .expect("-out opened at the object's start");
        let is_stdout = matches!(stream.rdbuf(), OutBuf::Cout(_));
        (stream, is_stdout)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
    /// The constructor of the object's formatter (`CBlastFormat`): the stream set to throw on
    /// `badbit`, which throws at once on a stream that failed before (`BLAST failed to write
    /// output`; the destructor of a formatter whose constructor threw does not run).
    fn begin_object(&mut self) -> anyhow::Result<()> {
        let (stream, _) = self.opened();
        stream.set_exceptions(IoState::BAD);
        if stream.take_thrown() {
            return Err(NativeError::write_failure().into());
        }
        self.formatting = true;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:318-324
    /// The formatter's destructor, at the object's end or on an error that leaves it: the
    /// stream's former mask (none) restored, then `flush`.
    fn end_object(&mut self) {
        if !std::mem::take(&mut self.formatting) {
            return;
        }
        let (stream, _) = self.opened();
        stream.set_exceptions(IoState::GOOD);
        stream.flush();
    }

    /// `cerr`'s tie to `cout` before a diagnostic line where the diagnostics go to `cerr`
    /// (`tied`): `cout` writes what it holds first (a file of `-out` does not). The posts of the
    /// loop come between objects, where the stream throws on nothing.
    fn before_diagnostic(&mut self, tied: bool) {
        if let Some(stream) = &mut self.stream {
            if tied && matches!(stream.rdbuf(), OutBuf::Cout(_)) {
                stream.flush();
                let _ = stream.take_thrown();
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:689-698
    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:310-323
    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:559-564
    /// Writes one serialized object (a Seq-annot or a Blast4-archive) into the stream through a
    /// fresh `ObjectOStream` (`writedb::metadata::ObjectOStream`, no diagnostic prefix): `write`
    /// receives the buffered sink and diagnostics stream, then `EndOfWrite`, and the destructor
    /// closes it. A close that fails posts `802.4` (`ObjectOStream::destroy`); a write that threw
    /// ends the run (`BLAST failed to write output`, exit status 6). Where the posts go to `cerr`
    /// and the stream is `cout`, the post's flush of `cout` (the tie) throws again inside the
    /// destructor: `std::terminate` (`NativeError::ios_failure_abort`), without the post.
    pub(super) fn write_object(
        &mut self,
        diag: &mut Diag,
        write: impl FnOnce(&mut dyn Write, &mut dyn Write) -> anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:6476-6492,7570-7610
        let tied = diag.is_cerr();
        let (stream, is_stdout) = self.opened();
        let tie = tied && is_stdout;
        let mut object = ObjectOStream::new(OStreamBuffer::new(stream), None);
        let mut sink = AnnotSink {
            buffer: &mut object.output,
            tie,
        };
        let written = write(&mut sink, &mut DiagWriter(&mut *diag)).and_then(|()| {
            object
                .output
                .end_of_write()
                .map_err(|_| anyhow::Error::from(thrown()))
        });
        let mut posts = Vec::new();
        object.destroy(&mut posts);
        if !posts.is_empty() {
            if tie {
                return Err(NativeError::ios_failure_abort().into());
            }
            diag.write_text(&posts);
        }
        written.map_err(|error| {
            if is_thrown(&error) {
                NativeError::write_failure().into()
            } else {
                error
            }
        })
    }

    /// The insertion of the Seq-annot of `desc` and `aligns` in `format` into the stream.
    fn write_annot(
        &mut self,
        format: SerialFormat,
        desc: AnnotDesc<'_>,
        aligns: &[V<'_>],
        buf: &mut Vec<u8>,
        diag: &mut Diag,
    ) -> anyhow::Result<()> {
        self.write_object(diag, |sink, diagnostics| {
            write_streamed(format, desc, aligns, buf, sink, diagnostics)
        })
    }
}

/// The `ios_base::failure` of the stream (a state change that met its exceptions mask), as the
/// I/O error of the writes of a Seq-annot (`AnnotSink`).
#[derive(Debug)]
struct Thrown;

impl std::fmt::Display for Thrown {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(crate::writedb::ofstream::IOS_FAILURE_WHAT)
    }
}

impl std::error::Error for Thrown {}

fn thrown() -> std::io::Error {
    std::io::Error::other(Thrown)
}

/// Whether `error` is the stream's throw (`thrown`).
fn is_thrown(error: &anyhow::Error) -> bool {
    error
        .downcast_ref::<std::io::Error>()
        .and_then(|error| error.get_ref())
        .is_some_and(|inner| inner.is::<Thrown>())
}

/// The bytes of a Seq-annot as its object stream's buffer takes them (`COStreamBuffer::Write`);
/// a throw of the stream is `thrown`. `flush` is `cerr`'s tie to `cout` before a post of a
/// replaced byte (`tie`): what the buffer and `cout` hold is written; otherwise nothing.
struct AnnotSink<'b, 'a, 'o> {
    buffer: &'b mut OStreamBuffer<'a, OutBuf<'o>>,
    tie: bool,
}

impl Write for AnnotSink<'_, '_, '_> {
    fn write(&mut self, data: &[u8]) -> std::io::Result<usize> {
        self.buffer.write_bytes(data).map_err(|_| thrown())?;
        Ok(data.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        if self.tie {
            self.buffer.end_of_write().map_err(|_| thrown())?;
            let stream = self.buffer.output();
            stream.flush();
            if stream.take_thrown() {
                return Err(thrown());
            }
        }
        Ok(())
    }
}

/// The diagnostics as a stream, for the posts of the bytes that NCBI's object stream replaces
/// while it writes a Seq-annot (`write_document_posting`).
struct DiagWriter<'d, 'a>(&'d mut Diag<'a>);

impl Write for DiagWriter<'_, '_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.write_text(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The database that an earlier object opened: its name, residue type and title, and the
/// warnings of its taxid list, which NCBI posts each time it opens it again.
struct OpenedDatabase {
    name: Vec<u8>,
    protein: bool,
    title: Vec<u8>,
    warnings: Vec<u8>,
}

impl<'r, 'o> Run<'r, 'o> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:796-818
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1657-1665
    /// The `CRemoteBlast` made of a stream and the loop of `LoadFromArchive`: the encoding guessed
    /// (`CFormatGuess`, accepting ASN.1 text, ASN.1 binary and ASN.1-style XML), then each object
    /// read, its messages posted and, unless it is the error archive, formatted.
    fn objects(&mut self, archive: &mut dyn Read, diag: &mut Diag) -> anyhow::Result<()> {
        let guessed = guess_encoding(archive);
        let mut reader = match guessed.encoding {
            Ok(Encoding::Text) => ArchiveReader::new(guessed.into_text(Names::BLAST_ARCHIVE)),
            Ok(Encoding::Binary) => {
                ArchiveReader::binary(guessed.into_binary(Names::BLAST_ARCHIVE))
            }
            Ok(Encoding::Xml) => {
                ArchiveReader::xml(guessed.into_xml(XmlStyle::Asn, Names::BLAST_ARCHIVE))
            }
            Ok(Encoding::Other) => return Err(invalid_archive()),
            Err(GuessStop::Exception(exception)) => return Err(exception_error(&exception)),
            Err(GuessStop::Undecided(what)) => return Err(Refusal::Undecided(what).into()),
        };
        while let Some((read, posts)) = reader.next_archive() {
            for post in &posts {
                post_line(&mut self.output, diag, &post.line(None));
            }
            let archive = read.map_err(read_error)?;
            self.object(archive, diag)?;
        }
        Ok(())
    }

    /// One object as `LoadFromArchive` reads it (its messages, then the request: `CImportStrategy`
    /// and, for the first search object, `GetSearchOptions`) and `PrintFormattedOutput` formats
    /// it.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1727-1739
    fn object(&mut self, archive: Archive, diag: &mut Diag) -> anyhow::Result<()> {
        self.post_messages(&archive.messages, diag)?;
        let search = match archive.search.map_err(stop_error)? {
            Some(search) => search,
            None => return Ok(()),
        };
        if let Subject::Database(name) = &search.request.subject {
            self.database = Some((name.clone(), search.request.db_is_protein()));
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1733-1737
        if search.request.service == b"psi" {
            self.psi_step = search
                .request
                .param(ParamSet::Format, b"Web_StepNumber")
                .map(|param| match param.value().choice() {
                    ("integer", value) => Ok(value.int() as i32),
                    _ => Err(()),
                });
        }
        if self.search_options.is_none() {
            let options = archive_request(&search.request).map_err(|e| exception_error(&e))?;
            self.hitlist = options.hitlist_size();
            self.algorithm_options = Some(RemoteOptions::built(&options));
            self.search_options = Some(options);
        }
        let printed = self.print(&search, diag);
        self.output.end_object();
        printed
    }

    /// The stored messages as `LoadFromArchive` posts them at the diagnostic level `Warning`:
    /// the text alone, in LOSAT's words (`wording::losat_bytes`); `Info` and `Trace` are below the
    /// level. A `Fatal` message ends NCBI's program (`abort`): refused.
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:505
    // (`LoadFromArchive` posts each message with `LOG_POST`, cited at
    // `report::archive::read::messages`.)
    fn post_messages(&mut self, messages: &[Message], diag: &mut Diag) -> anyhow::Result<()> {
        for message in messages {
            match message.severity {
                Severity::Info | Severity::Trace => {}
                Severity::Fatal => return Err(Refusal::FatalMessage.into()),
                Severity::Warning | Severity::Error | Severity::Critical => {
                    let mut line = crate::wording::losat_bytes(&message.text).into_owned();
                    line.push(b'\n');
                    post_line(&mut self.output, diag, &line);
                }
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:335-359
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:426-430
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:453-467
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2957-2962
    // (the formats up to 4: `-max_target_seqs`, or the larger of the description and alignment
    // counts where `-num_descriptions` or `-num_alignments` is given, cited at
    // `blastinput::report_counts::report_counts_from`.)
    /// Whether the command line's options set the hit list size of `ExtractAlgorithmOptions` for
    /// the output format `outfmt`, so that the options' own size is not kept.
    fn options_set_hitlist(&self, outfmt: i32) -> bool {
        let options = self.options;
        options.max_target_seqs.is_some()
            || options.num_alignments.is_some()
            || (outfmt <= output_format::FLAT_QUERY_ANCHORED_NO_IDENTITIES
                && options.num_descriptions.is_some())
    }

    /// `PrintFormattedOutput` of one search object, in NCBI's order: `-out` opened, the
    /// formatting options of the object (`extract_formatting` with the kept hit list size, its
    /// warnings), the check of SAM, the queries (`check_queries`), the database or the subject
    /// sequences (`subject_database`), the formatter's constructor (`Output::begin_object`), the
    /// check of the XInclude formats (13 and 14 without `-out`); the output formats that LOSAT
    /// does not write are refused there, after NCBI's errors that come before any output. The
    /// format 11 then writes the archive of the object with the hit list size set in its algorithm
    /// options (`CBlastFormat::WriteArchive`, `archive_out::write_archive`); the others read the
    /// results (`GetResultSet`: the alignments of each query, the masks), the prefetch's
    /// variable, and write the Seq-annot of each query that has alignments
    /// (`x_PrintStructuredReport`, `Output::write_annot`).
    fn print(&mut self, search: &Search, diag: &mut Diag) -> anyhow::Result<()> {
        self.output.open(&self.options.out)?;
        // A negative hit list size (an archive that NCBI's searches do not write, probes
        // P9.005, P9.009): `ExtractAlgorithmOptions` keeps it where no option sets the size, and
        // posts its warning (`hitlist_size < 5`) as for 0, which the handler reads in its place.
        let formatting = extract_formatting(self.options, self.hitlist.max(0) as usize)?;
        post_line(&mut self.output, diag, &formatting.warnings);
        let set_hitlist = self.options_set_hitlist(formatting.choice.number);
        if self.hitlist >= 0 || set_hitlist {
            self.hitlist = i32::try_from(formatting.counts.hitlist_size)
                .expect("a hit list size that the handler read from an int");
        }
        let search_options = self
            .search_options
            .as_ref()
            .expect("the options of the first search object");
        let program = program_type(search_options.program);
        sam_check(&formatting.choice, program.blastn)?;
        let request = &search.request;
        check_queries(&request.queries)?;
        let database = match &self.database {
            Some((name, protein)) => Some(opened_database_checked(
                &mut self.opened,
                name,
                *protein,
                search_options,
                self.settings,
            )),
            None => {
                check_subjects(&request.subject)?;
                None
            }
        };
        let database = match database {
            Some((outcome, warnings)) => {
                post_line(&mut self.output, diag, &warnings);
                Some(outcome?)
            }
            None => None,
        };
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:402-422
        self.output.begin_object()?;
        xinclude_check_for(&formatting.choice, Some(Path::new(&self.options.out)))?;
        if let Some(refusal) = refusals::output_format(formatting.choice.number) {
            return Err(refusal.into());
        }
        if formatting.choice.number == output_format::ARCHIVE_FORMAT {
            let algo = self
                .algorithm_options
                .as_mut()
                .expect("the options of the first search object");
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2978
            // `ExtractAlgorithmOptions` sets the hit list size in the options of every object
            // (`SetHitlistSize`), which only this format writes.
            algo.set(fields::HITLIST_SIZE, ParamValue::Integer(self.hitlist));
            let database = self.database.as_ref().map(|(name, _)| name.as_slice());
            return archive_out::write_archive(
                search,
                search_options,
                algo,
                database,
                self.psi_step,
                &mut self.output,
                &mut self.buf,
                diag,
            );
        }
        let format = match formatting.choice.number {
            output_format::ASN_TEXT => SerialFormat::AsnText,
            output_format::ASN_BINARY => SerialFormat::AsnBinary,
            _ => SerialFormat::Json,
        };
        let sets = search
            .results
            .alignment_sets(&request.queries)
            .map_err(stop_error)?;
        search
            .results
            .masks(request.queries.num_queries(), || {
                request.query_is_translated()
            })
            .map_err(stop_error)?;
        // `BlastFormatter_PreFetchSequenceData`: a value of `PRE_FETCH_SEQS_LIMIT` that is no
        // integer stops NCBI here (the search programs' one reader of the variable, decision
        // F3-6 (4)); NCBI's formats 8, 9 and 12 fetch nothing for any other value but `INT_MAX`
        // (`s_IsPrefetchFormat`, cited at `common::database::prefetch_reads_seq_aligns`), with
        // which NCBI loads each subject through its data loaders: refused (decision C3-1).
        if prefetch_seqs_limit()? == Some(i32::MAX) {
            return Err(Refusal::Prefetch.into());
        }
        let alignments = search.results.alignments();
        // `PruneSeqalign` takes the size as an unsigned int: a negative size is its value modulo
        // 2^32.
        let number = self.hitlist as u32;
        for set in sets {
            let set = &alignments[set];
            if set.is_empty() {
                continue;
            }
            let pruned = prune_seqalign(set, number).map_err(stop_error)?;
            let (task, program_number) = blast_type(search_options.program)?;
            let desc = AnnotDesc {
                task,
                program: program_number,
                database: match &database {
                    Some(database) => DatabaseDesc::of_name(
                        &database.name,
                        &database.title,
                        program.subject_nucleotide,
                    ),
                    None => DatabaseDesc::None,
                },
            };
            self.output
                .write_annot(format, desc, &pruned, &mut self.buf, diag)?;
        }
        Ok(())
    }
}

/// Posts `text` (whole diagnostic lines, none when empty) after the standard output written
/// before it (`Output::before_diagnostic`).
fn post_line(output: &mut Output<'_>, diag: &mut Diag, text: &[u8]) {
    if text.is_empty() {
        return;
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:6476-6492,7570-7610
    output.before_diagnostic(diag.is_cerr());
    diag.write_text(text);
}

/// The Seq-annot of one query of an archive: the program's type and database in its `desc`, and
/// the query's Seq-aligns as read, cut by `prune_seqalign`; `spill` takes the encoder's final
/// bytes after each Seq-align (`Spill`).
struct ArchiveDocument<'d> {
    desc: AnnotDesc<'d>,
    aligns: &'d [V<'d>],
    spill: Option<&'d RefCell<Spill<'d>>>,
}

impl SeqAnnotDocument for ArchiveDocument<'_> {
    fn desc(&self) -> AnnotDesc<'_> {
        self.desc
    }

    /// Each Seq-align as the archive holds it (`asn_types::write::write_value`): NCBI's
    /// formatter writes the objects it read.
    fn write_aligns<E: Encoder>(&self, enc: &mut E) {
        for align in self.aligns {
            enc.element();
            write_value(enc, align.ty, align.value);
            if let Some(spill) = self.spill {
                spill.borrow_mut().after_element(enc);
            }
        }
    }
}

/// The bytes that a Seq-annot's encoder has made final, handed on after a Seq-align once it holds
/// `FLUSH_BYTES` of them (`Encoder::flush_to`), so that a Seq-annot of any size is written with
/// little memory, as NCBI's object stream writes through its buffer. The ASN.1 text and binary
/// formats keep the last `SPILL_LAG` bytes back (`lag`): a Seq-annot that their plain encoder
/// stops on is written again through `Traced` behind the `written` bytes (`write_streamed`), and
/// NCBI's buffer (4096 bytes) had written more than those when it posted the first replaced byte.
pub(super) struct Spill<'w> {
    out: &'w mut dyn Write,
    lag: usize,
    /// The final bytes held back.
    held: Vec<u8>,
    /// The bytes handed to `out`.
    written: usize,
    /// The first error of `out`, after which nothing more is handed on.
    error: Option<std::io::Error>,
}

/// The bytes that `Spill` keeps back for the ASN.1 text and binary formats.
pub(super) const SPILL_LAG: usize = 2 * crate::writedb::ofstream::BUFSIZ;

impl<'w> Spill<'w> {
    pub(super) fn new(out: &'w mut dyn Write, lag: usize) -> Self {
        Spill {
            out,
            lag,
            held: Vec::new(),
            written: 0,
            error: None,
        }
    }

    pub(super) fn into_parts(self) -> (Vec<u8>, usize, Option<std::io::Error>) {
        (self.held, self.written, self.error)
    }

    /// After an element of a long list: the encoder's final bytes handed on, once it
    /// holds `FLUSH_BYTES` and has met no error.
    pub(super) fn after_element<E: Encoder>(&mut self, enc: &mut E) {
        if self.error.is_some() || enc.error().is_some() || enc.buffered() < FLUSH_BYTES {
            return;
        }
        let handed = if self.lag == 0 {
            let mut counted = Counted {
                inner: &mut *self.out,
                bytes: 0,
            };
            let handed = enc.flush_to(&mut counted, false);
            self.written += counted.bytes;
            handed
        } else {
            enc.flush_to(&mut self.held, false).and_then(|()| {
                let count = self.held.len().saturating_sub(self.lag);
                self.out.write_all(&self.held[..count])?;
                self.held.drain(..count);
                self.written += count;
                Ok(())
            })
        };
        if let Err(error) = handed {
            self.error = Some(error);
        }
    }
}

/// A writer that counts the bytes it passes on.
struct Counted<'w> {
    inner: &'w mut dyn Write,
    bytes: usize,
}

impl Write for Counted<'_> {
    fn write(&mut self, data: &[u8]) -> std::io::Result<usize> {
        let count = self.inner.write(data)?;
        self.bytes += count;
        Ok(count)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

/// A writer that drops the first `skip` bytes it is given (those that `Spill` handed on).
pub(super) struct SkipWriter<'w> {
    pub(super) inner: &'w mut dyn Write,
    pub(super) skip: usize,
}

impl Write for SkipWriter<'_> {
    fn write(&mut self, data: &[u8]) -> std::io::Result<usize> {
        let dropped = self.skip.min(data.len());
        self.skip -= dropped;
        self.inner.write_all(&data[dropped..])?;
        Ok(data.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

/// `write_document_to` of the Seq-annot of `desc` and `aligns` to `out`, element by element
/// (`Spill`): the plain encoder of `format`, then what it holds; a Seq-annot with a string byte
/// that NCBI's ASN.1 text or binary stream writes as `#` with an error message is written again
/// through `Traced` (`write_document_posting`, no prefix: blast_formatter sets none), the bytes
/// already handed on dropped. On another error of the encoder the bytes handed on stay written.
fn write_streamed(
    format: SerialFormat,
    desc: AnnotDesc<'_>,
    aligns: &[V<'_>],
    buf: &mut Vec<u8>,
    out: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> anyhow::Result<()> {
    buf.clear();
    let (plain, held, written, error) = {
        let spill = RefCell::new(Spill {
            out: &mut *out,
            lag: if format == SerialFormat::Json {
                0
            } else {
                SPILL_LAG
            },
            held: Vec::new(),
            written: 0,
            error: None,
        });
        let document = ArchiveDocument {
            desc,
            aligns,
            spill: Some(&spill),
        };
        let plain = write_document(format, &document, buf);
        let mut spill = spill.borrow_mut();
        (
            plain,
            std::mem::take(&mut spill.held),
            spill.written,
            spill.error.take(),
        )
    };
    if let Some(error) = error {
        return Err(error.into());
    }
    match plain {
        Ok(()) => {
            out.write_all(&held)?;
            out.write_all(buf)?;
            Ok(())
        }
        Err(SerialError::NonVisibleChar(byte)) => {
            let document = ArchiveDocument {
                desc,
                aligns,
                spill: None,
            };
            let mut skip = SkipWriter {
                inner: out,
                skip: written,
            };
            write_document_posting(format, byte, &document, None, buf, &mut skip, diagnostics)
        }
        Err(error) => Err(unwritable(error, format)),
    }
}

/// The refusal of a value of the Seq-annot that the plain encoder of `format` stopped on, other
/// than a byte that NCBI writes as `#` (`Refusal::Unwritable`); an error of another kind is the
/// writer's.
fn unwritable(error: SerialError, format: SerialFormat) -> anyhow::Error {
    let what = match error {
        SerialError::JsonDouble(_) => "a NaN or an infinite real number",
        SerialError::OctetString => "an OCTET STRING",
        SerialError::Unwritten(what) => what,
        other => return SeqAlignError::Serial(other).into(),
    };
    Refusal::Unwritable(what, format.outfmt()).into()
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1014-1040
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
/// `PruneSeqalign` of a Seq-align set to `number`, the hit list size (`m_HitlistSize`): the
/// Seq-aligns of the first `number` subjects, a new subject where the subject's id (row 1) does
/// not match the one before; a `disc` Seq-align counts as one more and is kept without the
/// check. NCBI writes the Seq-annot of a query with alignments also when nothing is left.
fn prune_seqalign<'a>(set: &[V<'a>], number: u32) -> Result<Vec<V<'a>>, Stop> {
    let mut kept = Vec::with_capacity(set.len());
    let mut previous: Option<&ParsedSeqId> = None;
    let mut count: u32 = 0;
    for &align in set {
        if align.field(member::seq_align::SEGS).choice().0 == "disc" {
            count = count.wrapping_add(1);
        } else {
            let subject = seq_align_id(align, 1)?;
            if previous.is_none_or(|previous| compare(subject, previous) != Sic::Yes) {
                count = count.wrapping_add(1);
            }
            if count > number {
                break;
            }
            previous = Some(subject);
        }
        kept.push(align);
    }
    Ok(kept)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:170-180
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:199-221
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:164-165
// NCBI reference (598d8ae6): c++/include/algo/blast/api/sseqloc.hpp:261-266
/// `x_ExtractQueries` and the scope of its first query (`GetScope` of the first): NCBI's errors of
/// a PSSM without a query or with a set of them (`std::runtime_error`, exit status 255) and of a
/// location other than `whole` or `int` (`LOSAT engine error`, exit status 3); a Bioseq without
/// an id (where NCBI reads the first id through a null pointer) and a request without a query
/// (the scope of the first of an empty vector) are refused.
fn check_queries(queries: &Queries) -> anyhow::Result<()> {
    let count = match queries {
        Queries::Pssm(_) => {
            let query = queries.pssm().and_then(|pssm| {
                pssm.field(member::pssmwithparameters::PSSM)
                    .get(member::pssm::QUERY)
            });
            let Some(query) = query else {
                return Err(NativeError::std_exception("PSSM has no query").into());
            };
            match query.choice() {
                ("seq", bioseq) if bioseq.field(member::bioseq::ID).item(0).is_none() => {
                    return Err(Refusal::Undecided(QUERY_WITHOUT_ID).into())
                }
                ("seq", _) => 1,
                _ => {
                    return Err(NativeError::std_exception(
                        "Cannot have multiple queries in a PSSM",
                    )
                    .into())
                }
            }
        }
        Queries::SeqLocList(_) => {
            let mut count = 0;
            for loc in queries.seq_locs() {
                if !matches!(loc.choice().0, "whole" | "int") {
                    return Err(engine_error(
                        "Only whole or int typed seq_loc is supported for CBlastQueryVector",
                    ));
                }
                count += 1;
            }
            count
        }
        Queries::BioseqSet(_) => {
            let bioseqs = queries.bioseqs();
            if bioseqs.iter().any(|bioseq| bioseq.ids().next().is_none()) {
                return Err(Refusal::Undecided(QUERY_WITHOUT_ID).into());
            }
            bioseqs.len()
        }
    };
    if count == 0 {
        return Err(Refusal::Undecided(
            "the scope of the first query of a request without a query",
        )
        .into());
    }
    Ok(())
}

/// What NCBI dereferences for a query Bioseq without an id.
const QUERY_WITHOUT_ID: &str = "the first id of a query that has no id";

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:234-246
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:390-392
/// `s_ConvertSubjects2TSeqLocVector` of the object's subject sequences (no database object came
/// before): a Bioseq without its length is NCBI's unset member (invalid input); one without an
/// id NCBI ranks (`*seqid` of a null `CRef`) and no subject at all (`front()` of an empty vector;
/// the archive reader makes a `seq-loc-list` subject an empty list) are refused.
fn check_subjects(subject: &Subject) -> anyhow::Result<()> {
    let bioseqs = subject.bioseqs();
    for bioseq in &bioseqs {
        if bioseq.best_id().is_none() {
            return Err(Refusal::Undecided("the best id of a subject that has none").into());
        }
        if bioseq.length().is_none() {
            return Err(invalid_archive());
        }
    }
    if bioseqs.is_empty() {
        return Err(Refusal::Undecided("the first subject of a request without one").into());
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:374-385
// NCBI reference (598d8ae6): c++/src/algo/blast/api/uniform_search.cpp:100-106
// NCBI reference (598d8ae6): c++/src/algo/blast/api/uniform_search.cpp:124-131
// The formatter gives the search database a positive list where the options' taxid list is not
// empty, and a negative list where their negative taxid list is not empty; the second list of a
// database throws (`CBlastException` `eInvalidArgument`, which `CATCH_ALL` reports as an engine
// error, exit status 3), before the database is opened.
/// `opened_database` after the formatter's id lists: both taxid lists of the options are NCBI's
/// engine error.
fn opened_database_checked<'d>(
    opened: &'d mut Option<OpenedDatabase>,
    name: &[u8],
    protein: bool,
    options: &BuiltOptions,
    settings: &ApplicationSettings,
) -> (anyhow::Result<&'d OpenedDatabase>, Vec<u8>) {
    let set = |list: &Option<Vec<i32>>| list.as_ref().is_some_and(|list| !list.is_empty());
    if set(&options.taxid_list) && set(&options.negative_taxid_list) {
        let exception = NcbiException::new(SET_NEGATIVE_GI_LIST, ID_LISTS);
        return (Err(NativeError::exception(&exception).into()), Vec::new());
    }
    opened_database(opened, name, protein, options, settings)
}

/// Where `CSearchDatabase::SetNegativeGiList` throws.
const SET_NEGATIVE_GI_LIST: ThrowSite = ThrowSite {
    file: "src/algo/blast/api/uniform_search.cpp",
    line: 127,
    kind: "CBlastException::eInvalidArgument",
    function: Some("ncbi::blast::CSearchDatabase::SetNegativeGiList"),
};

/// The message of `SET_NEGATIVE_GI_LIST`.
const ID_LISTS: &str = "Cannot have more than one type of id list filtering.";

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:364-387
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:285-297
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:380-398
/// The database of a database search as `PrintFormattedOutput` opens it, through the search
/// programs' open path (`common::database::open_search_database`, the list of
/// `SeqDb::restrict`, the alias files' counts, the formatter's reopening of each name and its
/// filtering algorithm): the database (its title for the Seq-annot) and the warnings of the
/// taxid list, which come before the outcome. NCBI opens the database with the taxid list of the
/// first object's options (`GetSeqDb`); an error of NCBI's database reader there makes NCBI ask
/// its servers for the database (`CRemoteBlastDbDataLoader`): refused (`Refusal::Database`,
/// decisions F3-5 and F3-6 (1) of line RG). The negative taxid list is built of the positive one
/// (an empty list where the request has only a negative one; both lists are
/// `opened_database_checked`'s error), so it restricts nothing. The
/// formatter's reopening of each name without the list (`GetBlastDbInfo`) and the filtering
/// algorithm (`Filtering algorithm ID N is not supported.`) end with `LOSAT database error`
/// (NCBI's `BLAST Database error`), exit status 2. An object whose database was opened for an
/// earlier one reuses it: the options are the first object's, so NCBI opens it again the same
/// way and posts the same warnings.
fn opened_database<'d>(
    opened: &'d mut Option<OpenedDatabase>,
    name: &[u8],
    protein: bool,
    options: &BuiltOptions,
    settings: &ApplicationSettings,
) -> (anyhow::Result<&'d OpenedDatabase>, Vec<u8>) {
    let reused = opened
        .as_ref()
        .is_some_and(|database| database.name == name && database.protein == protein);
    if !reused {
        let mut warnings = Vec::new();
        match open_database(name, protein, options, settings, &mut warnings) {
            Ok(title) => {
                *opened = Some(OpenedDatabase {
                    name: name.to_vec(),
                    protein,
                    title,
                    warnings,
                });
            }
            Err(error) => return (Err(error), warnings),
        }
    }
    let database = opened.as_ref().expect("opened above or before");
    (Ok(database), database.warnings.clone())
}

/// `opened_database` of a database not opened before: its title (`pairwise_database`), the
/// warnings of its taxid list appended to `warnings`.
fn open_database(
    name: &[u8],
    protein: bool,
    options: &BuiltOptions,
    settings: &ApplicationSettings,
    warnings: &mut Vec<u8>,
) -> anyhow::Result<Vec<u8>> {
    let text = std::str::from_utf8(name)
        .map_err(|_| Refusal::Database(String::from_utf8_lossy(name).into_owned()))?;
    let seqtype = if protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    };
    let mut db: SeqDb = open_search_database(name, seqtype, settings)
        .map_err(|error| remote_database(&error, text))?;
    if let Some(taxids) = options
        .taxid_list
        .as_ref()
        .filter(|taxids| !taxids.is_empty())
    {
        let list = DbList::Taxids(TaxidList::Positive(taxids.iter().copied().collect()));
        let mut posted = Vec::new();
        let restricted = db.restrict(&list, &mut posted);
        for warning in posted {
            warnings.extend(list_warning_line(warning, None));
        }
        restricted.map_err(|error| remote_database(&error, text))?;
    }
    db.check_alias_values()
        .map_err(|error| remote_database(&error, text))?;
    // `GetSearchOptions` keeps the filtering algorithm of the options only when it is not -1
    // (remote_blast.cpp:2056-2059).
    let filtering = options.db_filtering_algorithm_id.filter(|&id| id != -1);
    db.read_name_records(name, filtering.is_some())
        .map_err(|error| database_error(error, PROGRAM))?;
    if let Some(id) = filtering {
        db.name_mask_details(id)
            .map_err(|error| name_mask_error(&error, PROGRAM))?;
    }
    Ok(pairwise_database(&db).definition)
}

/// An error of NCBI's database reader where `GetSeqDb` opens the database `name`: a
/// `CSeqDBException` makes NCBI ask its servers (refused); another exception of NCBI's leaves
/// for `CATCH_ALL`; a database that LOSAT does not read is LOSAT's refusal.
fn remote_database(error: &SeqDbError, name: &str) -> anyhow::Error {
    match error.ncbi_exception() {
        Some(exception) if !exception.kind().starts_with("CSeqDBException::") => {
            NativeError::exception(exception).into()
        }
        Some(_) => Refusal::Database(name.to_string()).into(),
        None if error.ncbi_message().is_some() => Refusal::Database(name.to_string()).into(),
        None => database_error_of(error, PROGRAM),
    }
}

/// `Invalid input format for BLAST Archive.`: the `CInputException` that `Run` throws for an
/// input of another format and for a `CSerialException` of the loop (`BLAST query/options
/// error`, exit status 1).
pub(super) fn invalid_archive() -> anyhow::Error {
    options_error_bytes(b"Invalid input format for BLAST Archive.")
}

/// The archive reader's error, as the loop ends with it: a `CSerialException` is
/// `invalid_archive`; the end of the text inside the root object (`CEofException`, a
/// `CIOException` that `CATCH_ALL` ends without a post) keeps the exit status 0.
fn read_error(error: ReadError) -> anyhow::Error {
    match error {
        ReadError::Serial(_) => invalid_archive(),
        ReadError::Eof => NativeError::silent_end().into(),
        ReadError::Exception(exception) => exception_error(&exception),
        ReadError::Std(what) => NativeError::std_exception(what).into(),
        ReadError::OutOfMemory => NativeError::out_of_memory().into(),
        ReadError::Unsupported(text) => Refusal::UnreadArchive(text).into(),
    }
}

/// Where NCBI's code stops on the archive (`Stop`), as the run ends there: the serial exceptions
/// of a CHOICE or an unset member are `invalid_archive` (the loop's `catch`); an exception is
/// `exception_error`'s; a null pointer (an element stated as `NULL` in a list too), an `abort`
/// and a read of memory the input does not decide are refused.
pub(super) fn stop_error(stop: Stop) -> anyhow::Error {
    match stop {
        Stop::Exception(exception) => exception_error(&exception),
        Stop::InvalidSelection { .. } | Stop::Unassigned { .. } => invalid_archive(),
        Stop::NullPointer(what) | Stop::Abort(what) | Stop::Undefined(what) => {
            Refusal::Undecided(what).into()
        }
        Stop::NullElement => Refusal::Undecided("an element stated as NULL in a list").into(),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-228
/// An exception of NCBI's that leaves the loop, by its class: a serial exception is
/// `invalid_archive`; the others are `NativeError::exception`'s (a `CBlastException`,
/// `CRemoteBlastException` too, by its code and message; a `CSeqDBException`; `Error: ` with
/// `what()`). The archive reader and the request's options raise no `CInputException` or
/// `CArgException`.
fn exception_error(exception: &NcbiException) -> anyhow::Error {
    let kind = exception.kind();
    let class = kind.split("::").next().unwrap_or(kind);
    match class {
        "CSerialException" => invalid_archive(),
        class if crate::wording::SERIAL_EXCEPTION_CLASSES.contains(&class) => invalid_archive(),
        _ => NativeError::exception(exception).into(),
    }
}

#[cfg(test)]
#[path = "archive_tests.rs"]
mod archive_tests;
