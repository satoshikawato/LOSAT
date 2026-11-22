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
use super::ids::{fill_id_list, EmptyId, IdLists, IdSet, SeqIdSet, TextMatchIds};
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
        let exception = config_exception(
            273,
            INCONSISTENT,
            "ncbi::CWinMaskConfig::s_DetermineAppType",
            "one of '-mk_counts', '-convert' or '-ustat <stat_file>' must be specified",
        );
        return report(Failure::Ncbi(exception), diag);
    }
    let options = options(args, settings);
    match masking_run(args, &options, stdin, out, diag) {
        Ok(status) => status,
        Err(failure) => report(failure, diag),
    }
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
        ustat: args.value("ustat").expect("-ustat is given"),
        input: args.value("in").expect("-in has a default"),
        infmt: args.value("infmt").expect("-infmt has a default"),
        outfmt: args.value("outfmt").expect("-outfmt has a default"),
        parse_seqids: args.flag("parse_seqids"),
        params,
        pct,
        dust_level: dust.then(|| {
            args.integer("dust_level")
                .expect("-dust_level has a default")
        }),
        blastdb_registry: match &settings.blastdb {
            Ok(value) => Ok(value.as_deref()),
            Err(reason) => Err(reason.as_str()),
        },
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
    let input = if options.infmt == "blastdb" {
        None
    } else {
        Some(open_input_file("in", options.input, stdin).map_err(|_| {
            Failure::Ncbi(config_exception(
                381,
                "CWinMaskConfigException::can not open input stream",
                CONFIG,
                options.input,
            ))
        })?)
    };
    let mut source = work::reader(options, input)?;
    let mut output = args.output_file("out", out).map_err(Failure::Arg)?;
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

/// The id lists of `-ids` and `-exclude_ids` (`CIdSet_TextMatch`, or `CIdSet_SeqId` with
/// `-text_match false`), or the configuration's exception. `text_match` is `-text_match`, which
/// `-mk_counts` sets to true.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:374-376,405-445
/// An id of no words in a list matched as text is NCBI's crash (`ids::EmptyId`): refused. The
/// ids that `CSeq_id` cannot read are posted as they are read (`SeqIdSet::insert`).
fn id_lists(
    args: &Args,
    infmt: &str,
    text_match: bool,
    diag: &mut Diag,
) -> Result<IdLists, Failure> {
    let ids = args.value("ids").unwrap_or_default();
    let exclude_ids = args.value("exclude_ids").unwrap_or_default();
    if !ids.is_empty() && !exclude_ids.is_empty() {
        return Err(Failure::Ncbi(config_exception(
            412,
            INCONSISTENT,
            CONFIG,
            "only one of -ids or -exclude_ids can be specified",
        )));
    }
    let mut list = |name: &str, line: u32| -> Result<Option<IdSet>, Failure> {
        if name.is_empty() {
            return Ok(None);
        }
        if text_match {
            let mut set = TextMatchIds::default();
            fill_id_list(name, |id| set.insert(id)).map_err(|EmptyId { id }| {
                Failure::Refused(unsupported(&format!(
                    "an id list with an empty id ('{}': NCBI BLAST+ indexes its id sets at -1 and crashes)",
                    String::from_utf8_lossy(&id)
                )))
            })?;
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
        let Ok(()) = fill_id_list(name, |id| {
            set.insert(id, diag);
            Ok::<(), std::convert::Infallible>(())
        });
        Ok(Some(IdSet::SeqId(set)))
    };
    Ok(IdLists {
        ids: list(ids, 424)?,
        exclude_ids: list(exclude_ids, 439)?,
    })
}
