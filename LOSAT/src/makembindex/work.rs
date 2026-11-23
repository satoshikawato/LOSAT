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
use crate::ncbi_app::files::{open_text_file, Cin};
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
    let registry = match &settings.blastdb {
        Ok(value) => Ok(value.as_deref()),
        Err(reason) => Err(reason.as_str()),
    };
    match run_index(args, registry, stdin, out, diag) {
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
    if let Some(stat) = args.value("stat") {
        options.stat_file_name = (!stat.is_empty()).then(|| stat.to_string());
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
    base: &str,
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
        diag.cerr(format!("creating {name}...").as_bytes());
        let diag_cell = RefCell::new(&mut *diag);
        let made = make_index(
            input,
            &name,
            start,
            &mut stop,
            options,
            &mut |line| diag_cell.borrow_mut().post_bytes(line),
            &mut |text| diag_cell.borrow_mut().cerr(text),
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
        }
        num_seq = num_seq.wrapping_add(stop.wrapping_sub(start));
        if start == stop {
            diag.cerr(b"removed (empty)\n");
            return Ok((num_seq, num_vol));
        }
        num_vol += 1;
        diag.cerr(b"done\n");
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
fn fasta_file(name: &str, titles: bool) -> Result<FastaStream<'static>, Failure> {
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
    registry: Result<Option<&str>, &str>,
    stdin: &'a mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> Result<i32, Failure> {
    let old_style = args
        .boolean("old_style_index")
        .expect("-old_style_index has a default");
    let options = index_options(args, diag);
    let iformat = args
        .value("iformat")
        .unwrap_or(if old_style { "fasta" } else { "blastdb" });
    if !old_style && iformat == "fasta" {
        diag.err_post("new style index requires input format 'blastdb'");
        return Ok(1);
    }
    let mut seqstream: Box<dyn SequenceStream + 'a> = if iformat == "fasta" {
        if args.value("db_mask").is_some() {
            diag.err_post("-db_mask requires -iformat blastdb");
            return Ok(1);
        }
        match args.value("input") {
            Some(name) => Box::new(fasta_file(name, options.idmap)?),
            None => Box::new(FastaStream::new(Box::new(stdin), options.idmap)),
        }
    } else {
        let Some(input) = args.value("input") else {
            diag.err_post("input format 'blastdb' requires -input option");
            return Ok(1);
        };
        if args.flag("show_filters") {
            let db = BlastDbStream::open(input, Filter::None, registry, false)?;
            let descriptions = db.descriptions()?;
            let _ = out.write_all(format!("{descriptions}\n").as_bytes());
            let _ = out.flush();
            return Ok(0);
        }
        if !old_style {
            return new_style(args, input, registry, &options, diag);
        }
        let filter = match args.value("db_mask") {
            Some(name) => Filter::Name(name),
            None => Filter::None,
        };
        Box::new(BlastDbStream::open(input, filter, registry, options.idmap)?)
    };
    let base = if args.flag("show_filters") {
        String::new()
    } else {
        args.as_string("output").map_err(Failure::Arg)?.to_string()
    };
    make_volumes(seqstream.as_mut(), &base, &options, diag)?;
    Ok(0)
}

/// The new style: each volume of the database `input` indexed apart, with its super header.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:269-347
fn new_style(
    args: &Args,
    input: &str,
    registry: Result<Option<&str>, &str>,
    options: &Options,
    diag: &mut Diag,
) -> Result<i32, Failure> {
    if args.value("output").is_some() {
        warning(diag, "option 'output' is ignored for new style indices");
    }
    let db_vols = BlastDbStream::open(input, Filter::None, registry, false)?.volume_paths(input)?;
    let mask = args.value("db_mask");
    for dbv_name in &db_vols {
        let filter = match mask {
            Some(name) => Filter::Name(name),
            None => Filter::None,
        };
        let mut seqstream = BlastDbStream::open(dbv_name, filter, registry, options.idmap)?;
        let vol_num_seq = seqstream.num_oids();
        let (num_seq, num_vol) = make_volumes(&mut seqstream, dbv_name, options, diag)?;
        if num_seq != vol_num_seq {
            diag.err_post(&format!(
                "number of sequence reported by BLAST database volume ({vol_num_seq}) is not the \
                 same as in the index ({num_seq})"
            ));
            return Ok(1);
        }
        SuperHeader { num_seq, num_vol }.save(&format!("{dbv_name}.shd"));
    }
    Ok(0)
}
