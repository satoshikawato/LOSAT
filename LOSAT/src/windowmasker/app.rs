//! windowmasker's `Run` (`CWinMaskApplication::Run`) on the arguments that NCBI's reader read:
//! the mode (`CWinMaskConfig::s_DetermineAppType`), then the configuration in NCBI's order
//! (`CWinMaskConfig::CWinMaskConfig`: for the masking the input, its reader, the writer and the
//! id lists; for `-mk_counts` the id lists alone; for `-convert` nothing) and the work of the
//! mode (`counts::mk_counts`, `counts::convert`, `work::run`).

use std::io::Write;

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::{open_input_file, Cin};
use crate::ncbi_app::Diag;
use crate::utils::winmask::{MaskerParams, PctSettings};

use super::counts;
use super::ids::{fill_id_list, IdLists, IdSet, SeqIdSet, TextMatchIds};
use super::work::{
    self, config_exception, report, unsupported, Failure, Options, CONFIG, INCONSISTENT,
};

/// `CWinMaskApplication::Run` on `args`; the exit status.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:261-279
/// `stdin` is the standard input, which `-in -` reads.
pub(super) fn run(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    if args.flag("mk_counts") {
        let infmt = args.value("infmt").expect("-infmt has a default");
        return match id_lists(args, infmt, true, diag) {
            Ok(lists) => counts::mk_counts(args, settings, &lists, out, diag),
            Err(failure) => report(failure, diag),
        };
    }
    if args.flag("convert") {
        return counts::convert(args, out, diag);
    }
    if !args.given("ustat") {
        return report(no_mode(), diag);
    }
    let options = options(args, settings);
    match masking_run(args, &options, stdin, out, diag) {
        Ok(status) => status,
        Err(failure) => report(failure, diag),
    }
}

// No NCBI counterpart: LOSAT Web's check of the words (line NP WW3a, D-WW3-10, 2026-11-23; review
// of WW4c, F1; line NP WW4f, 2026-11-28).
/// `run` up to where it reads an input past its start, through `run`'s own steps, whose failure is
/// the check's, as `run` reports it (LOSAT Web's validation of the words, `super::check_in`): the
/// mode (`s_DetermineAppType`); for `-mk_counts` the id lists' exclusion (`exclusive_id_lists`)
/// and, without an id list to read, its refusal of `-infmt seqids` (`counts::check_mk_counts`) and
/// its opening of `-out`, its first file (`counts::counts_sink`: a name that does not open loses
/// the counts without an error, as NCBI's `ofstream`; LOSAT Web's refusal of an output that it
/// does not write); for `-convert` its exception for `-in -` (`counts::check_convert`); for the
/// masking, its first steps (`masking_opened`: `-in` opened and its reader made, then `-out`
/// opened). Where the run reads the contents of a file first (an id list; `-convert`'s counts),
/// the rest is the run's. 0 for an argv that passes them, as `run` goes on to read its files. The
/// check's files are LOSAT Web's check's (`web_adapter::check_tool_in`), which make no file.
pub(super) fn check(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    let checked = if args.flag("mk_counts") {
        exclusive_id_lists(args).and_then(|()| {
            let lists = ["ids", "exclude_ids"]
                .iter()
                .any(|name| !args.os_value(name).unwrap_or_default().is_empty());
            if lists {
                // `id_lists` reads the list before `mk_counts` refuses -infmt seqids.
                Ok(())
            } else {
                counts::check_mk_counts(args)
                    .and_then(|()| counts::counts_sink(args, out).map(drop))
            }
        })
    } else if args.flag("convert") {
        counts::check_convert(args)
    } else if args.given("ustat") {
        masking_opened(args, &options(args, settings), stdin, out).map(drop)
    } else {
        Err(no_mode())
    };
    match checked {
        Ok(()) => 0,
        Err(failure) => report(failure, diag),
    }
}

/// The error of an argv of no mode (`CWinMaskConfig::s_DetermineAppType`).
fn no_mode() -> Failure {
    Failure::Ncbi(config_exception(
        273,
        INCONSISTENT,
        "ncbi::CWinMaskConfig::s_DetermineAppType",
        "one of '-mk_counts', '-convert' or '-ustat <stat_file>' must be specified",
    ))
}

/// The masking settings of `CWinMaskConfig` (`eGenerateMasks` or `eGenerateMasksWithDuster`).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:322-368
/// `window_size` is a `Uint1` (`-window 300` is 44), the scores `Uint4`; their constraints keep
/// them positive (`-t_extend` may be 0, which is the file's value as when it is not given).
fn options<'a>(args: &'a Args, settings: &'a ApplicationSettings) -> Options<'a> {
    let unsigned = |name| args.integer(name).map_or(0, |value| value as u32);
    let params = MaskerParams {
        window_size: args.integer("window").map_or(0, |value| value as u8),
        textend: unsigned("t_extend"),
        cutoff_score: unsigned("t_thres"),
        max_score: unsigned("t_high"),
        min_score: unsigned("t_low"),
        set_max_score: unsigned("set_t_high"),
        set_min_score: unsigned("set_t_low"),
    };
    let pct = PctSettings {
        min_pct: args.double("t_low_pct"),
        extend_pct: args.double("t_extend_pct"),
        thres_pct: args.double("t_thres_pct"),
        max_pct: args.double("t_high_pct"),
    };
    let dust = args.boolean("dust").expect("-dust has a default");
    Options {
        ustat: args.os_value("ustat").expect("-ustat is given"),
        input: args.os_value("in").expect("-in has a default"),
        infmt: args.value("infmt").expect("-infmt has a default"),
        outfmt: args.value("outfmt").expect("-outfmt has a default"),
        parse_seqids: args.flag("parse_seqids"),
        params,
        pct,
        dust_level: dust.then(|| {
            args.integer("dust_level")
                .expect("-dust_level has a default")
        }),
        database: settings.database(),
    }
}

/// `CWinMaskConfig`'s constructor, then `work::run`; the exit status.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:320-327,375-401
/// The exceptions of the constructor leave it before its destructor could delete the writer:
/// a `maskinfo_*` writer made then writes nothing (unlike one of `Run`, `work::run`).
fn masking_run(
    args: &Args,
    options: &Options<'_>,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> Result<i32, Failure> {
    let (mut source, mut output) = masking_opened(args, options, stdin, out)?;
    let writer = work::writer(options.outfmt);
    let text_match = args
        .boolean("text_match")
        .expect("-text_match has a default");
    let lists = id_lists(args, options.infmt, text_match, diag)?;
    if lists.match_titles() {
        if let Some(source) = source.as_mut() {
            source.make_titles();
        }
    }
    work::run(options, source, writer, &lists, &mut output, diag)
}

/// The masking's source and output (`masking_opened`).
type Opened<'a, 'o> = (
    Option<Box<dyn work::SequenceSource + 'a>>,
    crate::ncbi_app::files::OutputFile<'o>,
);

/// The first steps of the masking (`masking_run`), which LOSAT Web's check takes too (`check`;
/// review of WW4c, F1; line NP WW4f): `-in` opened (none for `-infmt blastdb`) and its reader made
/// (`work::reader`), then `-out` opened (`AsOutputFile`, whose `CArgException` the application
/// reports as an argument error) and LOSAT Web's refusal of an output that it does not write
/// (`OpenFiles::create`).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:375-401
fn masking_opened<'a, 'o>(
    args: &Args,
    options: &Options<'_>,
    stdin: &'a mut Cin<'_>,
    out: &'o mut dyn Write,
) -> Result<Opened<'a, 'o>, Failure> {
    let input = if options.infmt == "blastdb" {
        None
    } else {
        Some(open_input_file("in", options.input, stdin).map_err(|_| {
            Failure::Ncbi(config_exception(
                381,
                "CWinMaskConfigException::can not open input stream",
                CONFIG,
                crate::utils::os_bytes::os_str_bytes(options.input),
            ))
        })?)
    };
    let source = work::reader(options, input)?;
    let output = args.output_file("out", out).map_err(Failure::Arg)?;
    if let Some(reason) = output.refusal() {
        return Err(Failure::Refused(unsupported(reason)));
    }
    Ok((source, output))
}

/// The exception of `-ids` with `-exclude_ids`, which `CWinMaskConfig` throws before it reads
/// either list (`id_lists`).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:412
fn exclusive_id_lists(args: &Args) -> Result<(), Failure> {
    let given = |name| !args.os_value(name).unwrap_or_default().is_empty();
    if given("ids") && given("exclude_ids") {
        return Err(Failure::Ncbi(config_exception(
            412,
            INCONSISTENT,
            CONFIG,
            "only one of -ids or -exclude_ids can be specified",
        )));
    }
    Ok(())
}

/// The id lists of `-ids` and `-exclude_ids` (`CIdSet_TextMatch`, or `CIdSet_SeqId` with
/// `-text_match false`), or the configuration's exception. `text_match` is `-text_match`, which
/// `-mk_counts` sets to true.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:374-376,405-445
/// An id of no words in a list matched as text (`ids::EmptyId`), where NCBI posts its error and
/// crashes, is posted as it is read and skipped (`PD-LOSAT-NCBI-DEFECTS` approved exception 22),
/// as are the ids that `CSeq_id` cannot read (`SeqIdSet::insert`).
fn id_lists(
    args: &Args,
    infmt: &str,
    text_match: bool,
    diag: &mut Diag,
) -> Result<IdLists, Failure> {
    exclusive_id_lists(args)?;
    let ids = args.os_value("ids").unwrap_or_default();
    let exclude_ids = args.os_value("exclude_ids").unwrap_or_default();
    let mut list = |name: &std::ffi::OsStr, line: u32| -> Result<Option<IdSet>, Failure> {
        if name.is_empty() {
            return Ok(None);
        }
        if text_match {
            let mut set = TextMatchIds::default();
            fill_id_list(name, |id| {
                if let Err(empty) = set.insert(id) {
                    diag.err_post_bytes(&empty.post());
                }
            });
            return Ok(Some(IdSet::Text(set)));
        }
        if infmt != "blastdb" {
            return Err(Failure::Ncbi(config_exception(
                line,
                INCONSISTENT,
                CONFIG,
                "-text_match false can be used only with infmt blastdb",
            )));
        }
        let mut set = SeqIdSet::default();
        fill_id_list(name, |id| set.insert(id, diag));
        Ok(Some(IdSet::SeqId(set)))
    };
    Ok(IdLists {
        ids: list(ids, 424)?,
        exclude_ids: list(exclude_ids, 439)?,
    })
}
