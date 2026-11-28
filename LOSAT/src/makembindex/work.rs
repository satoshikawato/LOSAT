//! `CMkIndexApplication::Run`: the options of the index from the arguments, the sequence stream
//! (a FASTA file or the standard input, or a BLAST database), and the index volumes that
//! `crate::dbindex::factory::make_index` writes one after another until the input ends, each
//! announced on the standard error stream (`creating <name>...`, then `done`, or
//! `removed (empty)` for the call that found no sequence left).
//!
//! The new style (the default) indexes each volume of the database apart: its volumes are
//! `<volume path>.NN.idx`, then `<volume path>.shd` holds the volume's number of sequences and
//! of index volumes. The old style (`-old_style_index true`) writes `<-output>.NN.idx` of the
//! whole input and no super header.

use std::cell::RefCell;
use std::io::Write;

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::dbindex::factory::{make_index, Stop};
use crate::dbindex::format::{volume_name, SuperHeader};
use crate::dbindex::stream::{BlastDbStream, FastaStream, Filter, SequenceStream, StreamError};
use crate::dbindex::{unsupported, Options, SeqNum, MAX_STRIDE};
use crate::ncbi_app::args::{ArgException, Args};
use crate::ncbi_app::files::{open_text_file, Cin, Ostream};
use crate::ncbi_app::Diag;
use crate::ncbi_exception::{post_line, NcbiException, Severity, ThrowSite};

/// Why `Run` ended without its own exit status.
enum Failure {
    /// An exception that leaves `Run`: the application reports it (`Application's execution
    /// failed`, error 106.16), exit status 3.
    Ncbi(NcbiException),
    /// A `CArgException` thrown in `Run` (`AsString` of `-output` without a value).
    Arg(ArgException),
    /// LOSAT's refusal, the whole line after `Error: `.
    Refused(String),
    /// A write to `-show_filters`, `-stat`, or an index output file met a closed pipe where
    /// SIGPIPE ends NCBI's process (`ncbi_app::files::Ostream`;
    /// PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 10): end at once with exit status 0 and no
    /// message.
    ClosedPipe,
}

impl From<StreamError> for Failure {
    fn from(error: StreamError) -> Self {
        match error {
            StreamError::Ncbi(exception) => Failure::Ncbi(exception),
            StreamError::Refused(reason) => Failure::Refused(reason),
        }
    }
}

/// makembindex's `Run` with the reports of the application around it; the exit status.
pub(super) fn run(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    match run_index(args, settings.database(), stdin, out, diag) {
        Ok(status) => status,
        Err(Failure::Ncbi(exception)) => {
            diag.write_text(&exception.report_default_ex(
                106,
                16,
                "Application's execution failed",
            ));
            3
        }
        Err(Failure::Arg(error)) => diag.run_arg_exception(error),
        Err(Failure::Refused(reason)) => diag.refuse(&reason),
        Err(Failure::ClosedPipe) => 0,
    }
}

/// A post of severity Warning (makembindex posts from that severity on).
fn warning(diag: &mut Diag, message: &str) {
    diag.write_text(&post_line(
        Severity::Warning,
        None,
        None,
        message.as_bytes(),
    ));
}

/// The options of the index from the arguments, with their warnings.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:153-206
fn index_options(args: &Args, diag: &mut Diag) -> Options {
    let mut options = Options::default_options();
    // `-verbosity` sets `report_level`, which nothing reads.
    let positive = |name| args.integer(name).map(|value| value as u64);
    if let Some(volsize) = positive("volsize") {
        options.max_index_size = volsize;
    }
    if let Some(stat) = args.bytes("stat") {
        options.stat_file_name = (!stat.is_empty()).then(|| stat.into_owned());
    }
    if let Some(nmer) = positive("nmer") {
        options.hkey_width = nmer;
    }
    options.legacy = args.boolean("legacy").expect("-legacy has a default");
    options.idmap = args.boolean("idmap").expect("-idmap has a default");
    if let Some(stride) = positive("stride") {
        if options.legacy {
            warning(diag, "-stride has no effect upon legacy index creation");
        } else {
            options.stride = stride;
        }
    }
    if let Some(requested) = positive("ws_hint") {
        if options.legacy {
            warning(diag, "-ws_hint has no effect upon legacy index creation");
        } else {
            let minimum = options.hkey_width + options.stride - 1;
            options.ws_hint = if requested < minimum {
                warning(
                    diag,
                    &format!(
                        "-ws_hint requested is too low. Setting to the minimum value of {minimum}"
                    ),
                );
                minimum
            } else {
                requested
            };
        }
    }
    options
}

/// The volumes of `input` named `<base>.NN.idx` from volume 0 on, until a call finds no
/// sequence left; the number of sequences indexed and of volumes written.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:310-332
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:350-368
fn make_volumes(
    input: &mut dyn SequenceStream,
    base: &[u8],
    options: &Options,
    diag: &mut Diag,
) -> Result<(u32, u32), Failure> {
    let (mut num_seq, mut num_vol) = (0u32, 0u32);
    let mut stop: SeqNum = 0;
    let mut vol_num = 0u32;
    loop {
        let start = stop;
        stop = SeqNum::MAX;
        let name = volume_name(base, vol_num);
        vol_num += 1;
        diag.cerr(&[b"creating ".as_slice(), &name, b"..."].concat());
        // A text that met a closed pipe of the standard error or of `-logfile` ended NCBI's
        // process there (`Diag::closed`): before the volume's file is made.
        if diag.closed() {
            return Err(Failure::ClosedPipe);
        }
        let diag_cell = RefCell::new(&mut *diag);
        let made = make_index(
            input,
            &name,
            start,
            &mut stop,
            options,
            &mut |line| {
                let mut diag = diag_cell.borrow_mut();
                diag.post_bytes(line);
                !diag.closed()
            },
            &mut |text| {
                let mut diag = diag_cell.borrow_mut();
                diag.cerr(text);
                !diag.closed()
            },
        );
        match made {
            Ok(()) => {}
            Err(Stop::Stream(StreamError::Ncbi(exception))) => {
                return Err(Failure::Ncbi(exception))
            }
            Err(Stop::Stream(StreamError::Refused(reason))) => {
                diag.cerr(b"\n");
                return Err(Failure::Refused(reason));
            }
            Err(Stop::NoDefinedIndex) => {
                diag.cerr(b"\n");
                return Err(Failure::Refused(no_defined_index(options)));
            }
            Err(Stop::ClosedPipe) => return Err(Failure::ClosedPipe),
        }
        num_seq = num_seq.wrapping_add(stop.wrapping_sub(start));
        let empty = start == stop;
        if empty {
            diag.cerr(b"removed (empty)\n");
        } else {
            num_vol += 1;
            diag.cerr(b"done\n");
        }
        if diag.closed() {
            return Err(Failure::ClosedPipe);
        }
        if empty {
            return Ok((num_seq, num_vol));
        }
    }
}

/// LOSAT's refusal of a volume for which NCBI's makembindex gives no defined index: a stride
/// of 32768 or more (NCBI's offsets then come from a shift past an `int`'s width, undefined in
/// C++; its x86 build finishes with what that build computes), or a stride above `ws_hint -
/// nmer + 1` (a zero `mult_`: NCBI's `Save` does not end).
fn no_defined_index(options: &Options) -> String {
    if options.stride >= MAX_STRIDE {
        unsupported(&format!(
            "-stride {} (32768 or more, where NCBI's index has no defined result: its offsets \
             come from shifting an int past its width)",
            options.stride
        ))
    } else {
        unsupported(&format!(
            "-stride {} above -ws_hint {} minus -nmer {} plus 1 (NCBI's makembindex does not \
             finish such an index volume)",
            options.stride, options.ws_hint, options.hkey_width
        ))
    }
}

/// `CSequenceIStreamFasta` of `name`: the file, or `failed to open input stream`.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_fasta.cpp:55-71
fn fasta_file(name: &std::ffi::OsStr, titles: bool) -> Result<FastaStream<'static>, Failure> {
    match open_text_file(name) {
        Ok(file) => Ok(FastaStream::new(file, titles)),
        Err(_) => Err(Failure::Ncbi(NcbiException::new(
            ThrowSite {
                file: "src/algo/blast/dbindex/sequence_istream_fasta.cpp",
                line: 64,
                kind: "CSequenceIStream_Exception::I/O error",
                function: None,
            },
            "failed to open input stream",
        ))),
    }
}

/// `Run` itself; `Ok` holds its exit status (an `exit` of 1 after an error post included).
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:151-378
fn run_index<'a>(
    args: &Args,
    database: crate::blastinput::ncbi_environment::DatabaseSettings<'_>,
    stdin: &'a mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> Result<i32, Failure> {
    let old_style = args
        .boolean("old_style_index")
        .expect("-old_style_index has a default");
    let options = index_options(args, diag);
    // A warning that met a closed pipe ended NCBI's process there (`Diag::closed`).
    if diag.closed() {
        return Err(Failure::ClosedPipe);
    }
    let iformat = args
        .value("iformat")
        .unwrap_or(if old_style { "fasta" } else { "blastdb" });
    if !old_style && iformat == "fasta" {
        diag.err_post("new style index requires input format 'blastdb'");
        return Ok(1);
    }
    let mut seqstream: Box<dyn SequenceStream + 'a> = if iformat == "fasta" {
        if args.os_value("db_mask").is_some() {
            diag.err_post("-db_mask requires -iformat blastdb");
            return Ok(1);
        }
        match args.os_value("input") {
            Some(name) => Box::new(fasta_file(name, options.idmap)?),
            None => Box::new(FastaStream::new(Box::new(stdin), options.idmap)),
        }
    } else {
        let Some(input) = args.bytes("input") else {
            diag.err_post("input format 'blastdb' requires -input option");
            return Ok(1);
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:236-240
        if args.flag("show_filters") {
            let db = BlastDbStream::open(&input, Filter::None, database, false)?;
            let descriptions = db.descriptions()?;
            let mut out = Ostream::new(out);
            let _ = out.write_all(format!("{descriptions}\n").as_bytes());
            let _ = out.flush();
            if out.closed() {
                return Err(Failure::ClosedPipe);
            }
            return Ok(0);
        }
        if !old_style {
            return new_style(args, &input, database, &options, diag);
        }
        let mask = args.bytes("db_mask");
        let filter = match &mask {
            Some(name) => Filter::Name(name),
            None => Filter::None,
        };
        Box::new(BlastDbStream::open(
            &input,
            filter,
            database,
            options.idmap,
        )?)
    };
    let base = if args.flag("show_filters") {
        Vec::new()
    } else {
        args.bytes("output")
            .ok_or_else(|| Failure::Arg(crate::ncbi_app::args::no_value("output")))?
            .into_owned()
    };
    make_volumes(seqstream.as_mut(), &base, &options, diag)?;
    Ok(0)
}

/// The new style: each volume of the database `input` indexed apart, with its super header.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:269-347
fn new_style(
    args: &Args,
    input: &[u8],
    database: crate::blastinput::ncbi_environment::DatabaseSettings<'_>,
    options: &Options,
    diag: &mut Diag,
) -> Result<i32, Failure> {
    if args.os_value("output").is_some() {
        warning(diag, "option 'output' is ignored for new style indices");
        if diag.closed() {
            return Err(Failure::ClosedPipe);
        }
    }
    let db_vols = BlastDbStream::open(input, Filter::None, database, false)?.volume_paths(input)?;
    let mask = args.bytes("db_mask");
    for dbv_name in &db_vols {
        let filter = match &mask {
            Some(name) => Filter::Name(name),
            None => Filter::None,
        };
        let mut seqstream = BlastDbStream::open(dbv_name, filter, database, options.idmap)?;
        let vol_num_seq = seqstream.num_oids();
        let (num_seq, num_vol) = make_volumes(&mut seqstream, dbv_name, options, diag)?;
        if num_seq != vol_num_seq {
            diag.err_post(&format!(
                "number of sequence reported by BLAST database volume ({vol_num_seq}) is not the \
                 same as in the index ({num_seq})"
            ));
            return Ok(1);
        }
        let super_header = SuperHeader { num_seq, num_vol };
        if super_header.save(&[dbv_name.as_slice(), b".shd"].concat()) {
            return Err(Failure::ClosedPipe);
        }
    }
    Ok(0)
}
