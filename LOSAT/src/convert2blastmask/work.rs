//! `CConvert2BlastMaskApplication::Run` on the arguments that NCBI's reader read: the reader of
//! `-in` (`CMaskFromFasta`), the writer of `-outfmt` (`CMaskWriterBlastDbMaskInfo`) and the
//! masks of each record.

use std::io::{Read, Write};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::{ArgException, Args};
use crate::ncbi_app::files::{file_exception, open_input_file, Cin};
use crate::ncbi_app::Diag;
use crate::ncbi_exception::{err_post_x, NcbiException};
use crate::seqdb::mask::FilterProgram;
use crate::seqid::model::{ObjectId, ParsedSeqId};
use crate::seqid::rank::{best_rank_score, find_best};
use crate::seqmasks_io::run::{finish, serial_write_failure, write_failure, RunFailure};
use crate::seqmasks_io::writer::{MaskInfoWriter, MaskWriter, Serial, WriteError};
use crate::writedb::fasta_source::reader::{
    out_of_range_what, pack_overflow_what, FastaReader, ReadError, ReaderFlags,
};

/// How a run stops early.
enum Failure {
    /// An NCBI exception that `Run` catches: `Run` writes `e.what()` and a line end to `cerr`, exit
    /// status 1.
    Caught(NcbiException),
    /// A `std::exception` that leaves `Run`: an error posted (`ERR_POST_X`, subcode 7) with
    /// "Application's execution failed: " and `e.what()`, exit status 3 (ncbiapp.cpp:675-679, cited
    /// at dustmasker's `report`).
    Std(String),
    /// A case LOSAT does not reproduce: a whole message.
    Refused(String),
    /// A write found the pipe closed where SIGPIPE ends NCBI's process (`OutputFile::closed`): the
    /// run ends at once, exit status 0, no message.
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
fn unsupported(what: &str) -> String {
    crate::seqmasks_io::run::unsupported(what, super::PROGRAM)
}

/// A `CArgException` of `AsInputFile` or `AsOutputFile`, which `Run` catches.
fn file_failure(error: ArgException) -> Failure {
    match file_exception(&error) {
        Some(exception) => Failure::Caught(exception),
        None => Failure::Refused(unsupported(&format!(
            "the argument error '{}' inside the run",
            error.message_text()
        ))),
    }
}

/// A failed read of `-in`: the `std::ios_base::failure` of NCBI's stream, which leaves `Run`.
fn read_failure(error: std::io::Error) -> Failure {
    crate::seqmasks_io::run::read_failure(error, super::PROGRAM)
}

/// `CConvert2BlastMaskApplication::Run` on `args`; the exit status.
///
/// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:219-236
/// The writer (`unique_ptr`) is destroyed as an exception leaves the `try` block, before the
/// `catch` writes: a `maskinfo_*` writer writes its object first (`run_masks`). `cerr` is tied
/// to `cout`, so the output is flushed before `e.what()`. `x_AddCmdOptions` and the usage
/// report write nothing.
pub(super) fn run(
    args: &Args,
    _settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    report(run_masks(args, stdin, out, diag), diag)
}

// No NCBI counterpart: LOSAT Web's check of the words (line NP WW4c; line NP WW4f, 2026-11-28).
/// `Run` up to where it reads the input past its start (LOSAT Web's check, `super::check_in`): its
/// first steps (`opened`: `-in` opened and its reader made, then `-out` opened), whose failure is
/// the check's, as `run` reports it; 0 where they pass, as the run then goes on to read the input.
/// The check's files are LOSAT Web's check's (`web_adapter::check_tool_in`), which make no file.
pub(super) fn check(
    args: &Args,
    _settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    report(opened(args, stdin, out).map(drop), diag)
}

/// The exit status of `run_masks`'s `result`, with what `Run`'s `catch` and the application write
/// for its failure (`run`).
fn report(result: Result<(), Failure>, diag: &mut Diag) -> i32 {
    match result {
        Ok(()) => 0,
        Err(Failure::Caught(exception)) => {
            diag.cerr(&[exception.what().as_slice(), b"\n"].concat());
            1
        }
        Err(Failure::Std(what)) => {
            diag.write_text(
                err_post_x(106, 7, &format!("Application's execution failed: {what}")).as_bytes(),
            );
            3
        }
        Err(Failure::Refused(reason)) => diag.refuse(&reason),
        Err(Failure::ClosedPipe) => 0,
    }
}

/// The body of `Run`'s `try` block, and the writer's destruction after it.
fn run_masks(
    args: &Args,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> Result<(), Failure> {
    let bytes = |name| {
        args.bytes(name)
            .expect("a key with a default or a mandatory key")
    };
    let outfmt = args.value("outfmt").expect("-outfmt has a default");
    let (mut reader, mut out) = opened(args, stdin, out)?;
    let mut writer = writer(
        outfmt,
        &bytes("masking_algorithm"),
        bytes("masking_options").into_owned(),
    );
    // A failed write leaves NCBI's stream bad, and the `interval` lines go on silently (exit
    // 0): `out` drops the bytes after the first failed write.
    let result = (|| {
        while let Some(record) = reader.next_sequence()? {
            if !record.mask.is_empty() {
                writer.print_id(&record.best_id, &record.mask, &mut out)?;
                serial_write_failure::<Failure>(&out, outfmt, super::PROGRAM)?;
            }
        }
        Ok(())
    })();
    finish(
        MaskWriter::MaskInfo(writer),
        &mut out,
        diag,
        result,
        outfmt,
        super::PROGRAM,
    )
}

/// `x_GetWriter`: makes a `CMaskWriterBlastDbMaskInfo` from `output`, `format`, a 0, the program
/// `prog` and `masking_options`; the program is named by `-masking_algorithm` in lower case.
///
/// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:199-217
/// (`AsOutputFile` opens `-out` first: `run_masks` opens it before this. `NStr::ToLower` is
/// `tolower` of each byte, the C locale's.)
fn writer(outfmt: &str, algorithm: &[u8], options: Vec<u8>) -> MaskInfoWriter {
    let program = match algorithm.to_ascii_lowercase().as_slice() {
        b"not_set" => FilterProgram::NotSet,
        b"dust" => FilterProgram::Dust,
        b"seg" => FilterProgram::Seg,
        b"windowmasker" => FilterProgram::Windowmasker,
        b"repeat" => FilterProgram::Repeat,
        _ => FilterProgram::Other,
    };
    let serial = match outfmt {
        "interval" => None,
        "maskinfo_asn1_bin" => Some(Serial::AsnBinary),
        "maskinfo_asn1_text" => Some(Serial::AsnText),
        "maskinfo_xml" => Some(Serial::Xml),
        other => unreachable!("-outfmt {other} passed its constraint"),
    };
    MaskInfoWriter::new(serial, 0, program, options)
}

/// One record of `CMaskFromFasta`: `GetBestID()` and `GetMask()`.
struct MaskedRecord {
    best_id: ParsedSeqId,
    mask: Vec<(u32, u32)>,
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:192-217
/// The first steps of `run_masks`, which LOSAT Web's check takes too (`check`; review of WW4c, F1;
/// line NP WW4f): `-in` opened and its reader made (`MaskFromFasta::new`), then `-out` opened
/// (`AsOutputFile`, whose `CArgException` `Run` catches: `file_failure`), and LOSAT Web's refusal
/// of an output that it does not write (`OpenFiles::create`), as soon as `-out` is opened
/// (windowmasker's place and words).
fn opened<'a, 'o>(
    args: &Args,
    stdin: &'a mut Cin<'_>,
    out: &'o mut dyn Write,
) -> Result<(MaskFromFasta<'a>, crate::ncbi_app::files::OutputFile<'o>), Failure> {
    let input = args.os_value("in").expect("-in has a default");
    let reader = MaskFromFasta::new(
        Box::new(open_input_file("in", input, stdin).map_err(file_failure)?),
        args.flag("parse_seqids"),
    )?;
    let out = args.output_file("out", out).map_err(file_failure)?;
    if let Some(reason) = out.refusal() {
        return Err(Failure::Refused(unsupported(reason)));
    }
    Ok((reader, out))
}

/// `CMaskFromFasta`: NCBI's `CFastaReader` with its flags
/// (`ReaderFlags::mask_from_fasta`), which keeps the intervals of each record.
struct MaskFromFasta<'a> {
    /// `None` once the reading has ended (at the end of the input or by an exception).
    reader: Option<FastaReader<Box<dyn Read + 'a>>>,
}

impl<'a> MaskFromFasta<'a> {
    /// `x_GetReader`: a new `CMaskFromFasta` over the input file of `-in` (`AsInputFile`) and the
    /// `-parse_seqids` flag.
    /// The constructor of `CFastaReader` over a stream reads the line reader's first buffer
    /// (cited at `seqmasks_io::fasta_reader::MaskFastaReader::new`): a failed read throws
    /// there, before `-out` is opened.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:192-197
    fn new(input: Box<dyn Read + 'a>, parse_seqids: bool) -> Result<Self, Failure> {
        let mut reader = FastaReader::new(input, ReaderFlags::mask_from_fasta(parse_seqids));
        match reader.take_read_error() {
            Some(error) => Err(read_failure(error)),
            None => Ok(Self {
                reader: Some(reader),
            }),
        }
    }

    /// `GetNextSequence`, and `GetBestID` and `GetMask` of the record it read; `None` at the
    /// end of the input. Every exception of `ReadOneSeq` leaves it, `eEOF` too (input of blank
    /// and comment lines only: `CFastaReader: Expected defline`), as no `ReadSet` catches it.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:67-74
    /// `HasMask()` is whether `GetMask()` has an interval (each interval sets it). `GetBestID`
    /// is the id that `ParseDefLine` chose (`FindBestChoice` over `GetIDs()` with
    /// `CSeq_id::BestRank`, fasta.cpp:593): with `-parse_seqids` among the record's ids, else the
    /// local id that the reader generated.
    fn next_sequence(&mut self) -> Result<Option<MaskedRecord>, Failure> {
        let Some(reader) = &mut self.reader else {
            return Ok(None);
        };
        if reader.at_eof() {
            let read_error = reader.take_read_error();
            self.reader = None;
            return match read_error {
                Some(error) => Err(read_failure(error)),
                None => Ok(None),
            };
        }
        // convert2blastmask posts at NCBI's release level (Error): the reader's warnings are
        // not written.
        let result = reader.read_one_seq(&mut |_warning| {});
        if let Some(error) = reader.take_read_error() {
            self.reader = None;
            return Err(read_failure(error));
        }
        let parses_ids = reader.parses_ids();
        let mask = reader.take_masks();
        let record = result.map_err(|error| {
            self.reader = None;
            match error {
                ReadError::Parse(error) => Failure::Caught(error.exception()),
                ReadError::Ncbi(exception) => Failure::Caught(*exception),
                ReadError::OutOfMemory => Failure::Std("std::bad_alloc".to_string()),
                ReadError::OutOfRange { pos, size } => Failure::Std(out_of_range_what(pos, size)),
                ReadError::PackOverflow => Failure::Refused(unsupported(&pack_overflow_what())),
            }
        })?;
        let best_id = if parses_ids {
            find_best(&record.ids, best_rank_score)
                .expect("a record has an id")
                .clone()
        } else {
            ParsedSeqId::Local(ObjectId::Id(record.local_id as i32))
        };
        Ok(Some(MaskedRecord { best_id, mask }))
    }
}
