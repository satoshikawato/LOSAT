//! `CConvert2BlastMaskApplication::Run` on the arguments that NCBI's reader read: the reader of
//! `-in` (`CMaskFromFasta`), the writer of `-outfmt` (`CMaskWriterBlastDbMaskInfo`) and the
//! masks of each record.

use std::io::{Read, Write};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::{ArgException, Args};
use crate::ncbi_app::files::{file_exception, open_input_file, Cin, Ostream};
use crate::ncbi_app::Diag;
use crate::ncbi_exception::{err_post_x, NcbiException};
use crate::seqdb::mask::FilterProgram;
use crate::seqid::model::{ObjectId, ParsedSeqId};
use crate::seqid::rank::{best_rank_score, find_best};
use crate::seqmasks_io::run::{finish, serial_write_failure, write_failure, RunFailure};
use crate::seqmasks_io::writer::{MaskInfoWriter, MaskWriter, Serial, WriteError};
use crate::writedb::fasta_source::reader::{FastaReader, ReadError, ReaderFlags};

/// How a run stops early.
enum Failure {
    /// An NCBI exception that `Run` catches: `cerr << e.what() << endl`, exit status 1.
    Caught(NcbiException),
    /// A `std::exception` that leaves `Run`: `ERR_POST_X(7, "Application's execution failed: "
    /// << e.what())`, exit status 3 (ncbiapp.cpp:675-679, quoted at dustmasker's `report`).
    Std(String),
    /// A case LOSAT does not reproduce: a whole message.
    Refused(String),
    /// A write found the pipe closed where SIGPIPE ends NCBI's process (`Ostream::closed`): the
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

/// A `CArgException` of `AsInputFile` or `AsOutputFile`, which `Run` catches.
fn file_failure(error: ArgException) -> Failure {
    match file_exception(&error) {
        Some(exception) => Failure::Caught(exception),
        None => Failure::Refused(crate::seqmasks_io::run::unsupported(
            &format!("the argument error '{}' inside the run", error.message()),
            super::PROGRAM,
        )),
    }
}

/// A failed read of `-in`: the `std::ios_base::failure` of NCBI's stream, which leaves `Run`.
fn read_failure(error: std::io::Error) -> Failure {
    crate::seqmasks_io::run::read_failure(error, super::PROGRAM)
}

/// `CConvert2BlastMaskApplication::Run` on `args`; the exit status.
///
/// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:219-236
/// ```c++
/// int CConvert2BlastMaskApplication::Run(void) {
///     int retval = 0;
///
///     try {
///         unique_ptr<CMaskFromFasta> reader(x_GetReader());
///         unique_ptr<CMaskWriterBlastDbMaskInfo> writer(x_GetWriter());
///
///         while (reader->GetNextSequence()) {
///             if(reader->HasMask()) writer->Print(reader->GetBestID(), reader->GetMask());
///         }
///     } catch (const CException& e) {
///         cerr << e.what() << endl;
///         retval = 1;
///     }
///     x_AddCmdOptions();
///     m_UsageReport.AddParam(CBlastUsageReport::eExitStatus, retval);
///     return retval;
/// }
/// ```
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
    let mut posts = Vec::new();
    let result = run_masks(args, stdin, out, &mut posts);
    diag.write_text(&posts);
    match result {
        Ok(()) => 0,
        Err(Failure::Caught(exception)) => {
            diag.cerr(format!("{}\n", exception.what()).as_bytes());
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
/// `posts` gets the messages that the object's stream posts (`MaskWriter::finish`).
fn run_masks(
    args: &Args,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    posts: &mut Vec<u8>,
) -> Result<(), Failure> {
    let value = |name| {
        args.value(name)
            .expect("a key with a default or a mandatory key")
    };
    let outfmt = value("outfmt");
    let mut reader = MaskFromFasta::new(
        Box::new(open_input_file("in", value("in"), stdin).map_err(file_failure)?),
        args.flag("parse_seqids"),
    )?;
    let mut output = args.output_file("out", out).map_err(file_failure)?;
    let mut writer = writer(outfmt, value("masking_algorithm"), value("masking_options"));
    // A failed write leaves NCBI's stream bad, and the `interval` lines go on silently (exit
    // 0): `out` drops the bytes after the first failed write.
    let mut out = Ostream::new(&mut output);
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
        result,
        outfmt,
        super::PROGRAM,
        posts,
    )
}

/// `x_GetWriter`: `CMaskWriterBlastDbMaskInfo(output, format, 0, prog, masking_options)`, the
/// program named by `-masking_algorithm` in lower case.
///
/// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:199-217
/// ```c++
/// CMaskWriterBlastDbMaskInfo*
/// CConvert2BlastMaskApplication::x_GetWriter() {
///     const CArgs& args = GetArgs();
///     const string& format(args["outfmt"].AsString());
///     CNcbiOstream& output = args["out"].AsOutputFile();
///
///     string algo=args["masking_algorithm"].AsString();
///     NStr::ToLower(algo);
///     EBlast_filter_program prog;
///     if     (algo == "not_set"     ) prog = eBlast_filter_program_not_set;
///     else if(algo == "dust"        ) prog = eBlast_filter_program_dust;
///     else if(algo == "seg"         ) prog = eBlast_filter_program_seg;
///     else if(algo == "windowmasker") prog = eBlast_filter_program_windowmasker;
///     else if(algo == "repeat"      ) prog = eBlast_filter_program_repeat;
///     else                            prog = eBlast_filter_program_other;
///
///     return(new CMaskWriterBlastDbMaskInfo(output, format, 0,
///                       prog, args["masking_options"].AsString()));
/// }
/// ```
/// (`AsOutputFile` opens `-out` first: `run_masks` opens it before this. `NStr::ToLower` is
/// `tolower` of each byte, the C locale's.)
fn writer(outfmt: &str, algorithm: &str, options: &str) -> MaskInfoWriter {
    let program = match algorithm.to_ascii_lowercase().as_str() {
        "not_set" => FilterProgram::NotSet,
        "dust" => FilterProgram::Dust,
        "seg" => FilterProgram::Seg,
        "windowmasker" => FilterProgram::Windowmasker,
        "repeat" => FilterProgram::Repeat,
        _ => FilterProgram::Other,
    };
    let serial = match outfmt {
        "interval" => None,
        "maskinfo_asn1_bin" => Some(Serial::AsnBinary),
        "maskinfo_asn1_text" => Some(Serial::AsnText),
        "maskinfo_xml" => Some(Serial::Xml),
        other => unreachable!("-outfmt {other} passed its constraint"),
    };
    MaskInfoWriter::new(serial, 0, program, options.as_bytes().to_vec())
}

/// One record of `CMaskFromFasta`: `GetBestID()` and `GetMask()`.
struct MaskedRecord {
    best_id: ParsedSeqId,
    mask: Vec<(u32, u32)>,
}

/// `CMaskFromFasta`: NCBI's `CFastaReader` with its flags
/// (`ReaderFlags::mask_from_fasta`), which keeps the intervals of each record.
struct MaskFromFasta<'a> {
    /// `None` once the reading has ended (at the end of the input or by an exception).
    reader: Option<FastaReader<Box<dyn Read + 'a>>>,
}

impl<'a> MaskFromFasta<'a> {
    /// `x_GetReader`: `new CMaskFromFasta(args["in"].AsInputFile(), args["parse_seqids"])`.
    /// The constructor of `CFastaReader` over a stream reads the line reader's first buffer
    /// (quoted at `seqmasks_io::fasta_reader::MaskFastaReader::new`): a failed read throws
    /// there, before `-out` is opened.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:192-197
    /// ```c++
    /// CMaskFromFasta*
    /// CConvert2BlastMaskApplication::x_GetReader() {
    ///     const CArgs& args = GetArgs();
    ///     CNcbiIstream& input = args["in"].AsInputFile();
    ///     return(new CMaskFromFasta(input, args["parse_seqids"]));
    /// }
    /// ```
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
    /// ```c++
    ///     bool GetNextSequence() {
    ///         m_hasMask = false;
    ///         m_mask.clear();
    ///         if (AtEOF()) return false;
    ///         SaveMask();
    ///         ReadOneSeq();
    ///         return true;
    ///     }
    /// ```
    /// `HasMask()` is whether `GetMask()` has an interval (each interval sets it). `GetBestID`
    /// is the id that `ParseDefLine` chose (`FindBestChoice(GetIDs(), CSeq_id::BestRank)`,
    /// fasta.cpp:593): with `-parse_seqids` among the record's ids, else the local id that the
    /// reader generated.
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
