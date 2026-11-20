//! windowmasker's `Run` (`CWinMaskApplication::Run`) on the arguments that NCBI's reader read:
//! the mode (`CWinMaskConfig::s_DetermineAppType`; LOSAT refuses the making and the conversion
//! of counts), then the configuration in NCBI's order (`CWinMaskConfig::CWinMaskConfig`: the
//! input, its reader, the writer, the id lists) and the masking (`work::run`).

use std::io::Write;

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::{open_input_file, Cin};
use crate::ncbi_app::Diag;
use crate::utils::winmask::{MaskerParams, PctSettings};

use super::ids::{fill_id_list, EmptyId, IdLists};
use super::work::{
    self, config_exception, report, unsupported, Failure, Options, CONFIG, INCONSISTENT,
};

/// `CWinMaskApplication::Run` on `args`; the exit status.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:261-279
/// ```c++
///     if(type == eAny){
///         if(args["mk_counts"])
///             type = eComputeCounts;
///         else if(args["convert"])
///             type = eConvertCounts;
///         else if(args["ustat"])
///             type = eGenerateMasksWithDuster;
///         else
///             NCBI_THROW( CWinMaskConfigException, eInconsistentOptions,
///                         "one of '-mk_counts', '-convert' or '-ustat <stat_file>' "
///                         "must be specified" );
///     }
///
///     if(type == eGenerateMasksWithDuster && !args["dust"].AsBoolean())
///         type = eGenerateMasks;
/// ```
/// `stdin` is the standard input, which `-in -` reads.
pub(super) fn run(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    if args.flag("mk_counts") {
        return diag.refuse(&unsupported(
            "-mk_counts (the making of a unit counts file)",
        ));
    }
    if args.flag("convert") {
        return diag.refuse(&unsupported(
            "-convert (the conversion of a unit counts file)",
        ));
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
/// ```c++
///       lstat_name( app_type >= eGenerateMasks ? args["ustat"].AsString() : "" ),
///       t_low_pct( app_type != eConvertCounts && args["t_low_pct"] ? args["t_low_pct"].AsDouble() : -1.0 ),
///       t_extend_pct( ... args["t_extend_pct"] ... ),
///       t_thres_pct( ... args["t_thres_pct"] ... ),
///       t_high_pct( ... args["t_high_pct"] ... ),
///       textend( app_type >= eGenerateMasks && args["t_extend"] ? args["t_extend"].AsInteger() : 0 ),
///       cutoff_score( app_type >= eGenerateMasks && args["t_thres"] ? args["t_thres"].AsInteger() : 0 ),
///       max_score( app_type != eConvertCounts && args["t_high"] ? args["t_high"].AsInteger() : 0 ),
///       min_score( app_type != eConvertCounts && args["t_low"] ? args["t_low"].AsInteger() : 0 ),
///       window_size( app_type >= eGenerateMasks && args["window"] ? args["window"].AsInteger() : 0 ),
/// ...
///       dust_level( app_type == eGenerateMasksWithDuster ? args["dust_level"].AsInteger() : 0 ),
/// ...
///         set_max_score = args["set_t_high"]  ? args["set_t_high"].AsInteger()
///                                             : 0;
///         set_min_score = args["set_t_low"]   ? args["set_t_low"].AsInteger()
///                                             : 0;
/// ```
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
        dust_level: dust.then(|| args.integer("dust_level").expect("-dust_level has a default")),
        blastdb_registry: match &settings.blastdb {
            Ok(value) => Ok(value.as_deref()),
            Err(reason) => Err(reason.as_str()),
        },
    }
}

/// `CWinMaskConfig`'s constructor, then `work::run`; the exit status.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:320-327,375-401
/// ```c++
///       is( app_type >= eGenerateMasks && args[kInputFormat].AsString() != "blastdb"
///                                      && determine_input ?
///           ( !(args[kInput].AsString() == "-")
///             ? new CNcbiIfstream( args[kInput].AsString().c_str() )
///             : static_cast<CNcbiIstream*>(&NcbiCin) ) : NULL ), reader( NULL ), writer( NULL ),
/// ...
///         if( is && !*is )
///         {
///             NCBI_THROW( CWinMaskConfigException,
///                         eInputOpenFail,
///                         args[kInput].AsString() );
///         }
///
///         if(determine_input && iformatstr != "seqids"){
///             ... (the reader, quoted at `work::reader`)
///         }
///
///         writer = x_GetWriter(args);
/// ```
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
    let source = work::reader(options, input)?;
    let mut output = args.output_file("out", out).map_err(Failure::Arg)?;
    let writer = work::writer(options.outfmt);
    let lists = id_lists(args, options.infmt)?;
    work::run(options, source, writer, &lists, &mut output, diag)
}

/// The id lists of `-ids` and `-exclude_ids` (`CIdSet_TextMatch`), or the configuration's
/// exception.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:405-445
/// ```c++
///     string ids_file_name( args["ids"].AsString() );
///     string exclude_ids_file_name( args["exclude_ids"].AsString() );
///
///     if(    !ids_file_name.empty()
///         && !exclude_ids_file_name.empty() )
///     {
///         NCBI_THROW( CWinMaskConfigException, eInconsistentOptions,
///                     "only one of -ids or -exclude_ids can be specified" );
///     }
///
///     if( !ids_file_name.empty() ) {
///         if( text_match ) {
///             ids = new CIdSet_TextMatch;
///         }else {
///             if( iformatstr == "blastdb" )
///                 ids = new CIdSet_SeqId;
///             else
///                 NCBI_THROW( CWinMaskConfigException, eInconsistentOptions,
///                         "-text_match false can be used only with "
///                         + string(kInputFormat) + " blastdb" );
///         }
///
///         FillIdList( ids_file_name, *ids );
///     }
///
///     if( !exclude_ids_file_name.empty() ) {
///         ... (the same, line 437)
///         FillIdList( exclude_ids_file_name, *exclude_ids );
///     }
/// ```
/// `-text_match false` with `-infmt blastdb` matches the lists as Seq-ids (`CIdSet_SeqId`),
/// whose messages for an id NCBI cannot parse name NCBI's build path: LOSAT refuses it. An id
/// of no words is NCBI's crash (`ids::EmptyId`): refused.
fn id_lists(args: &Args, infmt: &str) -> Result<IdLists, Failure> {
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
    let text_match = args.boolean("text_match").expect("-text_match has a default");
    let list = |name: &str, line: u32| -> Result<Option<_>, Failure> {
        if name.is_empty() {
            return Ok(None);
        }
        if !text_match {
            if infmt != "blastdb" {
                return Err(Failure::Ncbi(config_exception(
                    line,
                    INCONSISTENT,
                    CONFIG,
                    "-text_match false can be used only with infmt blastdb",
                )));
            }
            return Err(Failure::Refused(unsupported(
                "-text_match false (ids matched as Seq-ids, whose messages name NCBI's build)",
            )));
        }
        fill_id_list(name).map(Some).map_err(|EmptyId { id }| {
            Failure::Refused(unsupported(&format!(
                "an id list with an empty id ('{}': NCBI BLAST+ indexes its id sets at -1 and crashes)",
                String::from_utf8_lossy(&id)
            )))
        })
    };
    Ok(IdLists {
        ids: list(ids, 424)?,
        exclude_ids: list(exclude_ids, 439)?,
    })
}
