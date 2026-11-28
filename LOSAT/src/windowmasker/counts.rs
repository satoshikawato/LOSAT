//! windowmasker's `-mk_counts` (`CWinMaskCountsGenerator`: the unit counts of the input, file
//! after file) and `-convert` (`CWinMaskCountsConverter`: a counts file in another format). The
//! counts generator, the counts writers and the readers are `utils::winmask`'s; here are the
//! application's parts: the configuration, the input files and their readers, the output that
//! NCBI's writers open by name (silently, also when they cannot), and the messages.

use std::ffi::OsStr;
use std::fs::File;
use std::io::Write;

use crate::algorithm::common::database::TaxdbRun;
use crate::report::run_files::{Created, OpenFiles};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::{
    input_lines, standard_output_file, tool_stream, LineEnds, OutputFile,
};
use crate::ncbi_app::Diag;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::seqmasks_io::bdb_reader::{unset_seq_entry, SeqEntryAsked};
use crate::seqmasks_io::run::post_read_texts;
use crate::utils::os_bytes::path_from_bytes;
use crate::utils::winmask::{
    convert_counts, generate, CountsFormat, CountsWriter, FinalizeError, GenCountsSettings,
    GenError, GenomeSource, MaskerParams, PctSettings, WinMaskError, DEFAULT_TH,
};

use super::checkdup::check_duplicates;
use super::ids::IdLists;
use super::work::{
    database_source, fasta_source, open_counts, report, unsupported, winmask_failure, Failure,
    Record, SequenceSource,
};

/// Where the counts go (`CSeqMaskerOstat`'s stream).
enum Sink<'o> {
    /// NCBI's stream: the standard output (`-out -`) or a file opened by name
    /// (`OpenFiles::create`); after a failed write the rest is lost, as NCBI's `ofstream` loses it
    /// (`OutputFile`).
    Out(OutputFile<'o>),
    /// A name that did not open: NCBI's writer writes to its failed `ofstream`, which loses
    /// everything without a message.
    Lost,
}

/// The output of a counts writer of `format` made with `-out` `name`
/// (`CSeqMaskerOstatFactory::create`): the stream constructor for `-`, which the binary
/// format starts with an extra word on (NCBI's insertion, written before the next post,
/// `ncbi_app::Diag`); otherwise the name constructor, whose `ofstream` makes or empties the file at
/// once (an ascii writer given an empty name writes to the standard output), a file of the run's
/// files (`OpenFiles::create`; LOSAT Web's name of the standard output writes there, without the
/// stream's extra word). `stdout` is the run's standard output (`files::standard_output_file`).
/// `Err` with LOSAT's refusal of a name that LOSAT Web does not write.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_ascii.cpp:69-77
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_bin.cpp:55-67
fn open_sink<'o>(
    name: &OsStr,
    format: CountsFormat,
    stdout: &'o mut dyn Write,
    files: &OpenFiles,
) -> Result<Sink<'o>, Failure> {
    if name == "-" || (name.is_empty() && format == CountsFormat::Ascii) {
        let mut out = standard_output_file(stdout, files);
        if name == "-" && out.write_all(format.stream_prefix()).is_err() && out.closed() {
            return Err(Failure::ClosedPipe);
        }
        return Ok(Sink::Out(out));
    }
    match files.create(name) {
        Ok(Created::File(file)) => Ok(Sink::Out(OutputFile::Stream(tool_stream(
            Some(file),
            files,
        )))),
        Ok(Created::Stdout) => Ok(Sink::Out(standard_output_file(stdout, files))),
        Ok(Created::Refused(reason)) => Err(Failure::Refused(unsupported(&reason))),
        Err(_) => Ok(Sink::Lost),
    }
}

/// `finalize` of `writer` into `sink`, its message (`optimizing the data structure`) posted.
fn finalize(writer: CountsWriter, sink: &mut Sink<'_>, diag: &mut Diag) -> Result<(), Failure> {
    let mut lost = std::io::sink();
    let out: &mut dyn Write = match sink {
        Sink::Out(out) => out,
        Sink::Lost => &mut lost,
    };
    let mut closed = false;
    let result = writer.finalize(out, &mut |line| {
        diag.post(line);
        closed |= diag.closed();
    });
    if closed {
        return Err(Failure::ClosedPipe);
    }
    result.map_err(|error| match error {
        FinalizeError::Ncbi(exception) => Failure::Ncbi(exception),
        FinalizeError::Refused(what) => Failure::Refused(unsupported(&what)),
        FinalizeError::ClosedPipe => Failure::ClosedPipe,
    })
}

/// The file list of `-in` (`CWinMaskCountsGenerator::operator()`): the names between commas
/// (`NStr::Split` without flags: empty names are kept, an empty value has none), or with
/// `-fa_list` the lines of the file `-in` that are not empty, as the C++ library's `getline`
/// reads them from a text stream: lines end at `\n`, and a carriage return stays in the name
/// except where the library reads a text file's `\r\n` as `\n` (Windows); a file that does not
/// open or read gives no more names.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:228-242
pub(crate) fn file_list(input: &OsStr, fa_list: bool) -> Vec<Vec<u8>> {
    if !fa_list {
        if input.is_empty() {
            return Vec::new();
        }
        return crate::utils::os_bytes::os_str_bytes(input)
            .split(|&byte| byte == b',')
            .map(<[u8]>::to_vec)
            .collect();
    }
    let mut bytes = Vec::new();
    if let Ok(mut file) = File::open(input) {
        let _ = std::io::Read::read_to_end(&mut file, &mut bytes);
    }
    let text_mode = if cfg!(windows) {
        LineEnds::Windows
    } else {
        LineEnds::Unix
    };
    input_lines(&bytes, text_mode)
        .into_iter()
        .map(<[u8]>::to_vec)
        .collect()
}

/// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:636
/// The first read of a FASTA file that did not open: `CStreamReader` answers `eRW_Error` for a
/// stream that is not good, and `CBufferedLineReader` throws `CIOException::eRead`.
fn read_error() -> NcbiException {
    NcbiException::new(
        ThrowSite {
            file: "src/util/line_reader.cpp",
            line: 636,
            kind: "CIOException::eRead",
            function: Some("ncbi::CBufferedLineReader::x_ReadBuffer"),
        },
        "Read error",
    )
}

/// The input of `-mk_counts`: every file of the list, read with `CWinMaskUtil::CInputBioseq_CI`
/// (a FASTA reader without `-parse_seqids`, or a database reader), each sequence that
/// `CWinMaskUtil::consider` keeps.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:171-215
pub(super) struct GenomeInput<'r, 'd, 'w> {
    /// The names of the files, as their bytes (a `-fa_list` line is not text).
    files: Vec<Vec<u8>>,
    database: bool,
    lists: &'r IdLists,
    /// The settings of the database reader (`ApplicationSettings::database`).
    reader_settings: crate::blastinput::ncbi_environment::DatabaseSettings<'r>,
    diag: &'d mut Diag<'w>,
    /// What the last read posted (`GenomeSource::repeat_posts`).
    posts: Vec<u8>,
    /// The taxonomy database of the databases read so far: NCBI's one per process
    /// (`CSeqDBTaxInfo`), whose texts of opening come once, at the first lookup of the run's
    /// first read (`TaxdbRun::sequential_after`).
    taxonomy: Option<TaxdbRun>,
}

/// The function that `GenomeInput::each_record` gives each sequence, with the diagnostics.
pub(super) type RecordFn<'f> = dyn FnMut(&Record, &mut Diag) -> Result<(), Failure> + 'f;

/// Posts `line` (`LOG_POST`) after NCBI's writes of the output so far (NCBI's `cerr` is tied to
/// `cout`, `Diag`); a closed pipe ends the run there, as SIGPIPE ends NCBI's.
pub(super) fn post_after_output(diag: &mut Diag, line: &[u8]) -> Result<(), Failure> {
    diag.post_bytes(line);
    if diag.closed() {
        return Err(Failure::ClosedPipe);
    }
    Ok(())
}

impl GenomeInput<'_, '_, '_> {
    /// One read of every file: `f` gets each sequence that has letters and that the id lists
    /// keep, with the output and the diagnostics. The posts of the read are kept.
    pub(super) fn each_record(&mut self, f: &mut RecordFn<'_>) -> Result<(), Failure> {
        self.diag.start_record();
        let result = self.read_files(f);
        self.posts = self.diag.take_record();
        result
    }

    fn read_files(&mut self, f: &mut RecordFn<'_>) -> Result<(), Failure> {
        for name in &self.files {
            let mut source: Box<dyn SequenceSource> = if self.database {
                database_source(name, self.reader_settings, false)?
            } else {
                let file = File::open(path_from_bytes(name.as_slice()))
                    .map_err(|_| Failure::Ncbi(read_error()))?;
                fasta_source(file, false)?
            };
            if let Some(taxonomy) = &self.taxonomy {
                source.follow_taxonomy_of(taxonomy);
            }
            if self.lists.match_titles() {
                source.make_titles();
            }
            let database = self.database;
            let result = (|| -> Result<(), Failure> {
                while let Some(record) = source.next_record(self.diag)? {
                    post_read_texts::<Failure>(&record.texts, self.diag)?;
                    // NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:191-205
                    // `CInputBioseq_CI` asks each Seq-entry of a database for its Bioseq: the
                    // empty one of an OID without a Bioseq ends the run there (`-mk_counts` and
                    // `-checkdup`; the masking skips it, `work::mask_input`).
                    if database && record.iupac.is_empty() {
                        return Err(Failure::Ncbi(unset_seq_entry(SeqEntryAsked::Seq)));
                    }
                    if record.iupac.is_empty() || !self.lists.consider(&record) {
                        continue;
                    }
                    f(&record, self.diag)?;
                }
                Ok(())
            })();
            if let Some(taxonomy) = source.taxonomy() {
                self.taxonomy = Some(taxonomy.sequential_after());
            }
            result?;
        }
        Ok(())
    }
}

impl GenomeSource for GenomeInput<'_, '_, '_> {
    type Error = Failure;

    fn read(&mut self, f: &mut dyn FnMut(&[u8])) -> Result<(), Failure> {
        self.each_record(&mut |record, _| {
            f(&record.iupac);
            Ok(())
        })
    }

    fn repeat_posts(&mut self) -> Result<(), Failure> {
        let posts = std::mem::take(&mut self.posts);
        self.diag.write_text(&posts);
        self.posts = posts;
        if self.diag.closed() {
            return Err(Failure::ClosedPipe);
        }
        Ok(())
    }

    fn log(&mut self, line: &str) -> Result<(), Failure> {
        post_after_output(self.diag, line.as_bytes())
    }
}

/// `-mk_counts`: the configuration's id lists (`lists`, made before), then
/// `CWinMaskCountsGenerator`; the exit status.
///
/// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:146-198
pub(super) fn mk_counts(
    args: &Args,
    settings: &ApplicationSettings,
    lists: &IdLists,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    match mk_counts_run(args, settings, lists, out, diag) {
        Ok(status) => status,
        Err(failure) => report(failure, diag),
    }
}

/// The refusal of `-mk_counts` before it reads its input (LOSAT Web's validation too,
/// `app::check`): `-infmt seqids`, whose sequences NCBI fetches over the network.
pub(super) fn check_mk_counts(args: &Args) -> Result<(), Failure> {
    if args.value("infmt") == Some("seqids") {
        return Err(Failure::Refused(unsupported(
            "-mk_counts with -infmt seqids (NCBI BLAST+ fetches each sequence from NCBI over the network)",
        )));
    }
    Ok(())
}

/// The settings of the generator from `args`; `to_file` when the counts go to a named file,
/// whose constructor alone takes the `-t_*_pct` percentages.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:178-224
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:346-353
fn generator_settings(args: &Args, to_file: bool) -> GenCountsSettings {
    let unsigned = |name| args.integer(name).map_or(0, |value| value as u32);
    let mut th = DEFAULT_TH;
    if to_file {
        let pct = [
            args.double("t_low_pct"),
            args.double("t_extend_pct"),
            args.double("t_thres_pct"),
            args.double("t_high_pct"),
        ];
        for (slot, value) in th.iter_mut().zip(pct) {
            if let Some(value) = value {
                *slot = value;
            }
        }
    }
    GenCountsSettings {
        unit_size: args.integer("unit").map_or(0, |value| value as u8),
        genome_size: args
            .integer("genome_size")
            .map_or(0, |value| i64::from(value) as u64),
        mem: unsigned("mem"),
        min_count: unsigned("t_low"),
        max_count: unsigned("t_high"),
        th,
    }
}

fn mk_counts_run(
    args: &Args,
    settings: &ApplicationSettings,
    lists: &IdLists,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> Result<i32, Failure> {
    check_mk_counts(args)?;
    let infmt = args.value("infmt").expect("-infmt has a default");
    make_counts(args, settings, lists, infmt, out, diag)
}

/// `CWinMaskCountsGenerator` of the input `-in` in `-infmt` `infmt`; `stdout` is the standard
/// output. The diagnostics of a failure come after NCBI's writes of the output (`Diag`'s tie; a
/// closed pipe ends the run there, without them), and the stream ends after the run
/// (`OpenFiles::finish_streams`).
fn make_counts(
    args: &Args,
    settings: &ApplicationSettings,
    lists: &IdLists,
    infmt: &str,
    stdout: &mut dyn Write,
    diag: &mut Diag,
) -> Result<i32, Failure> {
    let input = args.os_value("in").expect("-in has a default");
    let files = file_list(input, args.boolean("fa_list").unwrap_or(false));
    let database = infmt == "blastdb";
    let output = args.os_value("out").expect("-out has a default");
    let format = counts_format(args);
    let mut sink = open_sink(output, format, stdout, args.files())?;
    let meta = args.bytes("meta").unwrap_or_default();
    let mut writer = CountsWriter::new(format, &meta);
    let mut source = GenomeInput {
        files,
        database,
        lists,
        reader_settings: settings.database(),
        diag: &mut *diag,
        posts: Vec::new(),
        taxonomy: None,
    };
    if args.boolean("checkdup").unwrap_or(false) {
        check_duplicates(&mut source)?;
    }
    let gen_settings = generator_settings(args, output != "-");
    generate(&gen_settings, &mut source, &mut writer).map_err(|error| match error {
        GenError::Source(failure) => failure,
        GenError::EmptyGenome(exception) => Failure::Ncbi(exception),
        GenError::Refused(what) => Failure::Refused(unsupported(&what)),
    })?;
    drop(source);
    finalize(writer, &mut sink, diag)?;
    Ok(0)
}

/// The counts format of `-sformat` and `-smem` (`CWinMaskConfig::SFormat`).
fn counts_format(args: &Args) -> CountsFormat {
    CountsFormat::from_sformat(
        args.value("sformat").expect("-sformat has a default"),
        args.integer("smem").expect("-smem has a default"),
    )
}

/// `CWinMaskCountsConverter`'s exception for `-in -` at `line` (the stream constructor's for
/// `-out -`, the name constructor's otherwise).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_counts_converter.cpp:47-82,140-148
fn input_name_exception(line: u32) -> NcbiException {
    NcbiException::new(
        ThrowSite {
            file: "src/algo/winmask/win_mask_counts_converter.cpp",
            line,
            kind: "Exception::argument error",
            function: Some("ncbi::CWinMaskCountsConverter::CWinMaskCountsConverter"),
        },
        "input file name must be non-empty",
    )
}

/// `-convert`: `CWinMaskCountsConverter` of `-in` to `-out` in `-sformat`; the exit status.
///
/// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:121-140
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_counts_converter.cpp:85-136
pub(super) fn convert(args: &Args, out: &mut dyn Write, diag: &mut Diag) -> i32 {
    match convert_run(args, out, diag) {
        Ok(status) => status,
        Err(failure) => report(failure, diag),
    }
}

fn convert_run(args: &Args, out: &mut dyn Write, diag: &mut Diag) -> Result<i32, Failure> {
    convert_file(args, out, diag)
}

/// An error of `-convert` of the counts file `file` (`-in <name>`): LOSAT's refusals in the words
/// of the conversion, else as [`winmask_failure`] gives it.
fn convert_failure(file: &str, error: WinMaskError) -> Failure {
    match error {
        WinMaskError::PctLine { .. } => {
            Failure::Refused(unsupported(&format!("-convert of {file}: {error};")))
        }
        WinMaskError::UnitSize { unit } => Failure::Refused(unsupported(&format!(
            "-convert of {file}, whose units of {unit} letters NCBI BLAST+ converts without an error while its reverse complement of a unit shifts the 32-bit unit by 32 bits or more (undefined behaviour in C++),"
        ))),
        error => winmask_failure(file, error),
    }
}

/// The exception of `-convert` before it reads its counts (LOSAT Web's validation too,
/// `app::check`): `-in -`, which the converter's constructors refuse.
pub(super) fn check_convert(args: &Args) -> Result<(), Failure> {
    if args.os_value("in") == Some(OsStr::new("-")) {
        let line = if args.os_value("out") == Some(OsStr::new("-")) {
            77
        } else {
            55
        };
        return Err(Failure::Ncbi(input_name_exception(line)));
    }
    Ok(())
}

/// `CWinMaskCountsConverter` of `-in` to `-out`; `stdout` is the standard output. The
/// diagnostics of a failure come after NCBI's writes of the output, as in [`make_counts`].
fn convert_file(args: &Args, stdout: &mut dyn Write, diag: &mut Diag) -> Result<i32, Failure> {
    check_convert(args)?;
    let input = args.os_value("in").expect("-in has a default");
    let output = args.os_value("out").expect("-out has a default");
    diag.post("reading counts...");
    let file = format!("-in {}", input.to_string_lossy());
    let istat = open_counts(
        input,
        &file,
        &MaskerParams::default(),
        &PctSettings::default(),
        diag,
    )?;
    let count_map = istat
        .count_map()
        .map_err(|error| convert_failure(&file, error))?;
    let meta = args.bytes("meta").unwrap_or_default();
    let metadata = if meta.is_empty() { istat.note() } else { &meta };
    let format = counts_format(args);
    let mut writer = CountsWriter::new(format, metadata);
    let mut sink = open_sink(output, format, stdout, args.files())?;
    // The messages come after NCBI's writes of the output so far (the binary format's extra
    // word); a closed pipe ends the run at the first, as SIGPIPE ends NCBI's.
    let mut closed = false;
    convert_counts(&istat, count_map, &mut writer, &mut |line| {
        closed = closed || post_after_output(diag, line.as_bytes()).is_err();
    })
    .map_err(|error| convert_failure(&file, error))?;
    if closed {
        return Err(Failure::ClosedPipe);
    }
    drop(istat);
    finalize(writer, &mut sink, diag)?;
    Ok(0)
}
