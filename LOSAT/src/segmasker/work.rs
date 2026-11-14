//! `SegMaskerApplication::Run` on the arguments that NCBI's reader read: the masker, the reader of
//! `-infmt`, the writer of `-outfmt`, and the masks of each protein sequence.

use std::io::Write;

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::{ArgException, Args};
use crate::ncbi_app::files::{file_exception, open_input_file, Cin, Ostream};
use crate::ncbi_app::Diag;
use crate::ncbi_exception::{err_post_x, NcbiException};
use crate::report::pairwise::cpp_default_double;
use crate::seqdb::mask::FilterProgram;
use crate::seqdb::SeqType;
use crate::seqid::ParsedSeqId;
use crate::seqmasks_io::bdb_reader::{iupacaa_letters, prints_title, BdbError, MaskBdbReader};
use crate::seqmasks_io::fasta_reader::{MaskFastaReader, MaskReadError};
use crate::seqmasks_io::run::{finish, serial_write_failure, write_failure, RunFailure};
use crate::seqmasks_io::writer::{MaskWriter, MaskedBioseq, Serial, WriteError};
use crate::utils::matrix::aa_char_to_ncbistdaa;
use crate::utils::seg::{SegMasker, SegParams};

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
            error.message()
        ))),
    }
}

/// An error of the database reader: NCBI's exception, which `Run` catches and writes with the
/// place of its throw (`e.what()`), or LOSAT's refusal where LOSAT does not know the place.
fn database_failure(error: BdbError) -> Failure {
    match error {
        BdbError::Ncbi(exception) if !exception.places_known() => {
            Failure::Refused(unsupported(&format!(
                "the report of the database error '{}', which names the place in NCBI's source \
                 that throws it,",
                exception.get_msg()
            )))
        }
        BdbError::Ncbi(exception) => Failure::Caught(exception),
        BdbError::Refused(reason) => Failure::Refused(reason),
        BdbError::Unsupported(what) => Failure::Refused(unsupported(&what)),
    }
}

/// `SegMaskerApplication::Run` on `args`; the exit status.
///
/// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:206-249
/// ```c++
/// int SegMaskerApplication::Run(void)
/// {
///     int retval = 0;
///     const CArgs& args = GetArgs();
///
///     try {
///
///         CRef<CObjectManager> objmgr(CObjectManager::GetInstance());
///
///         CSegMasker masker(args["window"].AsInteger(),
///                           args["locut"].AsDouble(),
///                           args["hicut"].AsDouble());
///
///         CRef<CSeq_entry> seq_entry;
///         unique_ptr<CMaskReader> reader(x_GetReader());
///         unique_ptr<CMaskWriter> writer(x_GetWriter());
///
///         while ( (seq_entry = reader->GetNextSequence()).NotEmpty() ) {
///
///         	// Allow skipping of oid
///         	if(seq_entry->Which() == CSeq_entry::e_not_set)
///         		continue;
///
///             CScope scope(*objmgr);
///             CSeq_entry_Handle seh = scope.AddTopLevelSeqEntry(*seq_entry);
///             CBioseq_Handle bioseq_handle = seh.GetSeq();
///             CSeqVector sequence_data =
///                 bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Ncbi);
///             unique_ptr<CSegMasker::TMaskList> masks(masker(sequence_data));
///             writer->Print(bioseq_handle, *masks, GetArgs()["parse_seqids"]);
///             // writer->Print(bioseq_handle, *masks);
///
///         }
///
///     } catch (const CException& e) {
///         cerr << e.what() << endl;
///         retval = 1;
///     }
///
///     return retval;
/// }
/// ```
/// The writer (`unique_ptr`) is destroyed as an exception leaves the `try` block, before the
/// `catch` writes: a `maskinfo_*` writer writes its object first (`run_masks`). `cerr` is tied
/// to `cout`, so the output is flushed before `e.what()`.
pub(super) fn run(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    match run_masks(args, settings, stdin, out, diag) {
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
fn run_masks(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> Result<(), Failure> {
    let value = |name| args.value(name).expect("a key with a default");
    // `CSegMasker masker(args["window"].AsInteger(), args["locut"].AsDouble(), ...)`.
    let double = |name| args.double(name).expect("an eDouble key with a default");
    let window = args.integer("window").expect("-window has a default");
    let (locut, hicut) = (double("locut"), double("hicut"));
    let outfmt = value("outfmt");
    let parse_seqids = args.flag("parse_seqids");
    let masker = seg_masker(window, locut, hicut);
    let registry = match &settings.blastdb {
        Ok(value) => Ok(value.as_deref()),
        Err(reason) => Err(reason.as_str()),
    };
    let mut source = reader(
        value("infmt"),
        value("in"),
        parse_seqids,
        prints_title(outfmt),
        registry,
        stdin,
    )?;
    let mut output = args.output_file("out", out).map_err(file_failure)?;
    let mut writer = writer(outfmt, window, locut, hicut);
    // A failed write leaves NCBI's stream bad, and the text writers go on silently (exit 0):
    // `out` drops the bytes after the first failed write.
    let mut out = Ostream::new(&mut output);
    let result = (|| {
        while let Some(record) = source.next_record()? {
            let masks = masker.seq_buffer_seg::<false>(&record.stdaa);
            let mask: Vec<(u32, u32)> = masks
                .iter()
                .map(|interval| (interval.start as u32, (interval.end - 1) as u32))
                .collect();
            let seq = MaskedBioseq {
                ids: &record.ids,
                title: record.title.as_deref(),
                taxname: record.taxname.as_ref(),
                iupac: &record.iupacaa,
                protein: true,
            };
            writer.print(&seq, &mask, parse_seqids, &mut out, diag)?;
            serial_write_failure::<Failure>(&out, outfmt, super::PROGRAM)?;
        }
        Ok(())
    })();
    finish(writer, &mut out, diag, result, outfmt, super::PROGRAM)
}

// NCBI reference (598d8ae6): c++/src/algo/segmask/segmask.cpp:40-52
// ```c++
// CSegMasker::CSegMasker(int window /* = kSegWindow */,
//                        double locut /* = kSegLocut */,
//                        double hicut /* = kSegHicut */)
// : m_SegParameters(SegParametersNewAa())
// {
//     ...
//     m_SegParameters->window = window;
//     m_SegParameters->locut = locut;
//     m_SegParameters->hicut = hicut;
// }
// ```
/// `CSegMasker(window, locut, hicut)`: `SegParametersNewAa` with the three values as given;
/// `SeqBufferSeg` checks them (`s_SegParametersCheck`: `SegParams::new`, a window not above 0
/// is 12) and keeps `overlaps` `FALSE` (`SegMasker::seq_buffer_seg::<false>`). The `maxbogus`
/// of 2 is not lowered to a window of 1 as NCBI's check does: a window of one letter has at most
/// one, so no window has more than either.
fn seg_masker(window: i32, locut: f64, hicut: f64) -> SegMasker {
    let window = usize::try_from(window).unwrap_or(0);
    SegMasker::with_params(&SegParams::new(window, locut, hicut))
}

// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:181-190
// ```c++
// string BuildAlgorithmParametersString(const CArgs& args)
// {
//     CNcbiOstrstream os;
//     if (args.Exist("locut") &&
//         args.Exist("hicut") &&
//         args.Exist("window")) {
//         // SEG
//         os << "window=" << args["window"].AsInteger() << "; "
//            << "locut=" << args["locut"].AsDouble() << "; "
//            << "hicut=" << args["hicut"].AsDouble();
// ```
/// `x_GetWriter`: the writer of `-outfmt` (whose constraint leaves no other value).
///
/// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:170-200
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
///     } else if (NStr::StartsWith(format, "maskinfo_asn1_bin")) {
///         CNcbiOstream& output = args[kOutput].AsOutputFile(CArgValue::fBinary);
///         retval =
///             new CMaskWriterBlastDbMaskInfo(output, format, 1,
///                                eBlast_filter_program_seg,
///                                BuildAlgorithmParametersString(args));
///     } else if (NStr::StartsWith(format, "maskinfo_")) {
///         ... (the same)
/// ```
/// (`AsOutputFile` opens `-out` in each branch: `run_masks` opens it before this.)
fn writer(outfmt: &str, window: i32, locut: f64, hicut: f64) -> MaskWriter {
    let mask_info = |serial| {
        let parameters = format!(
            "window={window}; locut={}; hicut={}",
            cpp_default_double(locut),
            cpp_default_double(hicut)
        );
        MaskWriter::mask_info(Some(serial), 1, FilterProgram::Seg, parameters.into_bytes())
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

/// One protein of the input (`CBioseq_Handle`): its Seq-ids, its title, its NCBIstdaa codes
/// (`eCoding_Ncbi`, what SEG reads) and its IUPACaa letters (`eCoding_Iupac`, what the FASTA
/// writer writes).
struct Record {
    ids: Vec<ParsedSeqId>,
    title: Option<Vec<u8>>,
    /// The name of a database sequence's source, from which its empty title is made
    /// (`seqmasks_io::bdb_reader::BdbRecord::taxname`); none for a FASTA record.
    taxname: Option<std::sync::Arc<[u8]>>,
    stdaa: Vec<u8>,
    iupacaa: Vec<u8>,
}

/// `CMaskReader::GetNextSequence` of the FASTA reader or the database reader.
enum Source<'a> {
    Fasta(Box<MaskFastaReader<Box<dyn std::io::Read + 'a>>>),
    Database(Box<MaskBdbReader>),
}

/// `x_GetReader`: `CMaskFastaReader` of `-in` (a protein reader), or `CMaskBDBReader` of the
/// protein database `-in`.
///
/// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:152-168
/// ```c++
///     if (format == "fasta") {
///         CNcbiIstream& input = args[kInput].AsInputFile();
///         retval = new CMaskFastaReader(input, false, args["parse_seqids"]);
///     } else if (format == "blastdb") {
///         retval = new CMaskBDBReader(args[kInput].AsString(), false);
///     } else {
///         _ASSERT("Unknown input format" == 0);
///     }
/// ```
/// (`-infmt` takes only the two values of its constraint.)
fn reader<'a>(
    infmt: &str,
    input: &str,
    parse_seqids: bool,
    title_printed: bool,
    registry: Result<Option<&str>, &str>,
    stdin: &'a mut Cin<'_>,
) -> Result<Source<'a>, Failure> {
    match infmt {
        "fasta" => {
            let input: Box<dyn std::io::Read + 'a> =
                Box::new(open_input_file("in", input, stdin).map_err(file_failure)?);
            MaskFastaReader::new(input, true, parse_seqids)
                .map(|reader| Source::Fasta(Box::new(reader)))
                .map_err(read_failure)
        }
        "blastdb" => MaskBdbReader::open(
            input,
            SeqType::Protein,
            registry,
            title_printed,
            super::PROGRAM,
        )
        .map(|reader| Source::Database(Box::new(reader)))
        .map_err(database_failure),
        other => unreachable!("-infmt {other} passed its constraint"),
    }
}

/// A failed read of `-in`: the `std::ios_base::failure` of NCBI's stream, which leaves `Run`.
fn read_failure(error: std::io::Error) -> Failure {
    crate::seqmasks_io::run::read_failure(error, super::PROGRAM)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1480-1503 (`s_CompOn`)
// ```c
// 	Int4 letter;
// 	char* seq,* seqmax;
// ...
// 	while (seq < seqmax) {
// 		letter = *seq++;
// 		if (!alphaflag[letter])
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1683-1685 (`s_ShiftWin1`)
// ```c
//     if ((Uint1) win->seq[length] == FENCE_SENTRY) {
//         return FALSE;
//     }
// ```
// (`FENCE_SENTRY` is 201, blast_util.h:364.) A database byte of 128 or more indexes NCBI's table
// before its start (undefined); makeblastdb writes 0 to 27, and 28 to 127 are bogus letters to
// SEG in both. Review of line RG item K2, M1.
impl Source<'_> {
    fn next_record(&mut self) -> Result<Option<Record>, Failure> {
        match self {
            Source::Fasta(reader) => {
                // segmasker posts at NCBI's release level (Error): the reader's warnings are not
                // written.
                let record = reader
                    .next(&mut |_warning| {})
                    .map_err(|error| match error {
                        // A `CObjReaderParseException`, which `Run` catches.
                        MaskReadError::Ncbi(exception) => Failure::Caught(exception),
                        MaskReadError::Read(error) => read_failure(error),
                        MaskReadError::OutOfMemory => Failure::Std("std::bad_alloc".to_string()),
                        MaskReadError::Refused(what) => Failure::Refused(unsupported(&what)),
                    })?;
                Ok(record.map(|record| Record {
                    ids: record.ids,
                    title: record.title,
                    taxname: None,
                    // IUPACaa letters: the 20 of SEG's alphabet have their NCBIstdaa codes, and
                    // every other letter is one SEG counts as bogus, as NCBIstdaa's `*` and `-`
                    // that the IUPACaa `X` stands for.
                    stdaa: record
                        .letters
                        .iter()
                        .map(|&letter| aa_char_to_ncbistdaa(letter))
                        .collect(),
                    iupacaa: record.letters,
                }))
            }
            Source::Database(reader) => {
                let Some(record) = reader.next_record().map_err(database_failure)? else {
                    return Ok(None);
                };
                if let Some(&code) = record.residues.iter().find(|&&code| code >= 128) {
                    return Err(Failure::Refused(unsupported(&format!(
                        "a database residue byte {code}, above 127 (NCBI's SEG reads it as a \
                         negative char, and stops its window at 201),"
                    ))));
                }
                Ok(Some(Record {
                    ids: record.ids,
                    title: Some(record.title),
                    taxname: record.taxname,
                    iupacaa: iupacaa_letters(&record.residues),
                    stdaa: record.residues,
                }))
            }
        }
    }
}
