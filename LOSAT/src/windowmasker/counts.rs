//! windowmasker's `-mk_counts` (`CWinMaskCountsGenerator`: the unit counts of the input, file
//! after file) and `-convert` (`CWinMaskCountsConverter`: a counts file in another format). The
//! counts generator, the counts writers and the readers are `utils::winmask`'s; here are the
//! application's parts: the configuration, the input files and their readers, the output that
//! NCBI's writers open by name (silently, also when they cannot), and the messages.

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::Ostream;
use crate::ncbi_app::Diag;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::utils::winmask::{
    convert_counts, generate, CountsFormat, CountsWriter, FinalizeError, GenCountsSettings,
    GenError, GenomeSource, MaskerParams, PctSettings, DEFAULT_TH,
};

use super::checkdup::check_duplicates;
use super::ids::IdLists;
use super::work::{
    considered, database_source, fasta_source, open_counts, report, unsupported, winmask_failure,
    Failure, Record, SequenceSource,
};

/// Where the counts go (`CSeqMaskerOstat`'s stream).
enum Sink<'a> {
    /// The standard output (`-out -`).
    Stdout(Ostream<'a>),
    /// A file opened by name; after a failed write the rest is lost, as NCBI's `ofstream`
    /// loses it.
    File(Option<BufWriter<File>>),
    /// A name that did not open: NCBI's writer writes to its failed `ofstream`, which loses
    /// everything without a message.
    Lost,
}

impl Write for Sink<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Sink::Stdout(out) => out.write(buf),
            Sink::File(file) => {
                if let Some(writer) = file {
                    if writer.write_all(buf).is_err() {
                        *file = None;
                    }
                }
                Ok(buf.len())
            }
            Sink::Lost => Ok(buf.len()),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Sink::Stdout(out) => out.flush(),
            Sink::File(file) => {
                if let Some(writer) = file {
                    if writer.flush().is_err() {
                        *file = None;
                    }
                }
                Ok(())
            }
            Sink::Lost => Ok(()),
        }
    }
}

/// The output of a counts writer of `format` made with `-out` `name`
/// (`CSeqMaskerOstatFactory::create`): the stream constructor for `-`, which the binary
/// format starts with an extra word on; otherwise the name constructor, whose `ofstream` makes
/// or empties the file at once (an ascii writer given an empty name writes to the standard
/// output).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_ascii.cpp:69-77
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_bin.cpp:55-67
fn open_sink<'a>(name: &str, format: CountsFormat, stdout: &'a mut dyn Write) -> Sink<'a> {
    if name == "-" || (name.is_empty() && format == CountsFormat::Ascii) {
        let mut out = Ostream::new(stdout);
        if name == "-" && !format.stream_prefix().is_empty() {
            let _ = out.write_all(format.stream_prefix());
            let _ = out.flush();
        }
        return Sink::Stdout(out);
    }
    match File::create(name) {
        Ok(file) => Sink::File(Some(BufWriter::with_capacity(1 << 16, file))),
        Err(_) => Sink::Lost,
    }
}

/// `finalize` of `writer` into `sink`, its message (`optimizing the data structure`) posted.
fn finalize(writer: CountsWriter, sink: &mut Sink<'_>, diag: &mut Diag) -> Result<(), Failure> {
    writer
        .finalize(sink, &mut |line| diag.post(line))
        .map_err(|error| match error {
            FinalizeError::Ncbi(exception) => Failure::Ncbi(exception),
            FinalizeError::Refused(what) => Failure::Refused(unsupported(&what)),
            FinalizeError::Io(_) => Failure::ClosedPipe,
        })
}

/// The file list of `-in` (`CWinMaskCountsGenerator::operator()`): the names between commas
/// (`NStr::Split` without flags: empty names are kept, an empty value has none), or with
/// `-fa_list` the lines of the file `-in` that are not empty (`getline`: a carriage return
/// stays in the name; a file that does not open or read gives no more names).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:228-242
fn file_list(input: &str, fa_list: bool) -> Vec<String> {
    if !fa_list {
        if input.is_empty() {
            return Vec::new();
        }
        return input.split(',').map(str::to_string).collect();
    }
    let mut bytes = Vec::new();
    if let Ok(mut file) = File::open(input) {
        let _ = std::io::Read::read_to_end(&mut file, &mut bytes);
    }
    bytes
        .split(|&b| b == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| String::from_utf8_lossy(line).into_owned())
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
    files: Vec<String>,
    database: bool,
    lists: &'r IdLists,
    registry: Result<Option<&'r str>, &'r str>,
    pub(super) diag: &'d mut Diag<'w>,
    /// What the last read posted (`GenomeSource::repeat_posts`).
    posts: Vec<u8>,
}

impl GenomeInput<'_, '_, '_> {
    /// One read of every file: `f` gets each sequence that has letters and that the id lists
    /// keep, with the diagnostics. The posts of the read are kept.
    pub(super) fn each_record(
        &mut self,
        f: &mut dyn FnMut(&Record, &mut Diag),
    ) -> Result<(), Failure> {
        self.diag.start_record();
        let result = self.read_files(f);
        self.posts = self.diag.take_record();
        result
    }

    fn read_files(&mut self, f: &mut dyn FnMut(&Record, &mut Diag)) -> Result<(), Failure> {
        for name in &self.files {
            let mut source: Box<dyn SequenceSource> = if self.database {
                database_source(name, self.registry, false)?
            } else {
                let file = File::open(name).map_err(|_| Failure::Ncbi(read_error()))?;
                fasta_source(file, false)?
            };
            while let Some(record) = source.next_record(self.diag)? {
                for text in &record.texts {
                    self.diag.err_post(text);
                }
                if record.iupac.is_empty() || !considered(self.lists, &record)? {
                    continue;
                }
                f(&record, self.diag);
            }
        }
        Ok(())
    }
}

impl GenomeSource for GenomeInput<'_, '_, '_> {
    type Error = Failure;

    fn read(&mut self, f: &mut dyn FnMut(&[u8])) -> Result<(), Failure> {
        self.each_record(&mut |record, _| f(&record.iupac))
    }

    fn repeat_posts(&mut self) {
        let posts = std::mem::take(&mut self.posts);
        self.diag.write_text(&posts);
        self.posts = posts;
    }

    fn log(&mut self, line: &str) {
        self.diag.post(line);
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
    let infmt = args.value("infmt").expect("-infmt has a default");
    if infmt == "seqids" {
        return Err(Failure::Refused(unsupported(
            "-mk_counts with -infmt seqids (NCBI BLAST+ fetches each sequence from NCBI over the network)",
        )));
    }
    let output = args.value("out").expect("-out has a default");
    let format = counts_format(args);
    let mut sink = open_sink(output, format, out);
    let input = args.value("in").expect("-in has a default");
    let mut writer = CountsWriter::new(format, args.value("meta").unwrap_or_default().as_bytes());
    let mut source = GenomeInput {
        files: file_list(input, args.boolean("fa_list").unwrap_or(false)),
        database: infmt == "blastdb",
        lists,
        registry: match &settings.blastdb {
            Ok(value) => Ok(value.as_deref()),
            Err(reason) => Err(reason.as_str()),
        },
        diag: &mut *diag,
        posts: Vec::new(),
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
    let input = args.value("in").expect("-in has a default");
    let output = args.value("out").expect("-out has a default");
    if input == "-" {
        let line = if output == "-" { 77 } else { 55 };
        return Err(Failure::Ncbi(input_name_exception(line)));
    }
    diag.post("reading counts...");
    let file = format!("-in {input}");
    let istat = open_counts(
        input,
        &file,
        &MaskerParams::default(),
        &PctSettings::default(),
        diag,
    )?;
    let count_map = istat
        .count_map()
        .map_err(|error| winmask_failure(&file, error))?;
    let meta = args.value("meta").unwrap_or_default();
    let metadata = if meta.is_empty() {
        istat.note()
    } else {
        meta.as_bytes()
    };
    let format = counts_format(args);
    let mut writer = CountsWriter::new(format, metadata);
    let mut sink = open_sink(output, format, out);
    convert_counts(&istat, count_map, &mut writer, &mut |line| diag.post(line))
        .map_err(|error| winmask_failure(&file, error))?;
    drop(istat);
    finalize(writer, &mut sink, diag)?;
    Ok(0)
}
