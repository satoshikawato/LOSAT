//! `losat windowmasker`: NCBI BLAST+ 2.17.0 windowmasker 1.0.0 (`CWinMaskApplication`), the
//! masker of over-represented and low-complexity regions of nucleotide sequences by the counts
//! of their units (WindowMasker, `utils::winmask`) and, with `-dust`, the symmetric DUST
//! algorithm. Its command line is NCBI's (`ncbi_app`), not clap's: `main` hands it the words
//! after `windowmasker`. Its three modes are the making of a counts file (`-mk_counts`, with
//! `-checkdup` its check for duplicated sequences), the conversion of a counts file to another
//! format (`-convert`) and the masking with a counts file (`-ustat`). The work of its `Run` is in
//! `app` (the configuration, `CWinMaskConfig`), `counts` and `checkdup` (the counts) and `work`
//! (the masking).

mod app;
mod checkdup;
mod counts;
mod ids;
mod sdust;
mod work;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::ncbi_app::args::{flag, key, Constraint, Dep, Desc, Kind};
use crate::ncbi_app::Tool;
use crate::ncbi_version::Component;

/// The name in LOSAT's refusals.
const PROGRAM: &str = "WINDOWMASKER";

// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:76,99-110
/// The text of `windowmasker -help` (NCBI 2.17.0's, `CArgDescriptions::PrintUsage` of the
/// descriptions in `DESCS` after `HideStdArgs`), frozen. Its part before `OPTIONAL ARGUMENTS` is
/// the usage that `-h` and the argument errors print.
const HELP: &str = include_str!("help.txt");

const INT_MAX: i64 = i32::MAX as i64;

// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:50-213
// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_cmdline_args.cpp:40-51
/// The arguments of `CWinMaskConfig::AddWinMaskArgs(arg_desc)` (`eAny`, the input determined)
/// in the order it describes them; NCBI's reader adds the help flags before them and `-logfile`
/// and `-conffile` after them (`ncbi_app::args`). The values that a `CArgAllow_Strings` allows
/// are listed in the order of their `std::set`.
const DESCS: &[Desc] = &[
    key("ustat", Kind::Text, None, None),
    key("in", Kind::InputFile, Some("-"), None),
    key("out", Kind::OutputFile, Some("-"), None),
    key("checkdup", Kind::Boolean, Some("false"), None),
    key("fa_list", Kind::Boolean, Some("false"), None),
    key(
        "mem",
        Kind::Integer,
        Some("1536"),
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key("meta", Kind::Text, Some(""), None),
    key(
        "unit",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, 16)),
    ),
    key("genome_size", Kind::Integer, None, None),
    key(
        "window",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key(
        "t_extend",
        Kind::Integer,
        None,
        Some(Constraint::Integers(0, INT_MAX)),
    ),
    key(
        "t_thres",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key(
        "set_t_high",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key(
        "set_t_low",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    flag("parse_seqids"),
    key(
        "outfmt",
        Kind::Text,
        Some("interval"),
        Some(Constraint::Strings(&[
            "fasta",
            "interval",
            "maskinfo_asn1_bin",
            "maskinfo_asn1_text",
            "maskinfo_xml",
            "seqloc_asn1_bin",
            "seqloc_asn1_text",
            "seqloc_xml",
        ])),
    ),
    key(
        "t_high",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key("t_high_pct", Kind::Double, None, Some(PCT)),
    key("t_thres_pct", Kind::Double, None, Some(PCT)),
    key("t_extend_pct", Kind::Double, None, Some(PCT)),
    key(
        "t_low",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key("t_low_pct", Kind::Double, None, Some(PCT)),
    key(
        "infmt",
        Kind::Text,
        Some("fasta"),
        Some(Constraint::Strings(&["blastdb", "fasta", "seqids"])),
    ),
    key("exclude_ids", Kind::Text, Some(""), None),
    key("ids", Kind::Text, Some(""), None),
    key("text_match", Kind::Boolean, Some("T"), None),
    key(
        "sformat",
        Kind::Text,
        Some("ascii"),
        Some(Constraint::Strings(&[
            "ascii", "binary", "oascii", "obinary",
        ])),
    ),
    key("smem", Kind::Integer, Some("512"), None),
    key("dust", Kind::Boolean, Some("F"), None),
    key("dust_level", Kind::Integer, Some("20"), None),
    flag("mk_counts"),
    flag("convert"),
];

/// `CArgAllow_Doubles( 0.0, 100.0 )` of the `-t_*_pct` keys.
const PCT: Constraint = Constraint::Doubles(0.0, 100.0, "0..100");

// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:212-254
/// The `SetDependency` calls of `AddWinMaskArgs`, in their order.
fn dependencies() -> Vec<(&'static str, Dep, &'static str)> {
    let excludes = |of: &'static str, names: &[&'static str]| {
        names
            .iter()
            .map(move |&name| (of, Dep::Excludes, name))
            .collect::<Vec<_>>()
    };
    let mut list = excludes(
        "mk_counts",
        &[
            "outfmt",
            "ustat",
            "window",
            "t_thres",
            "t_extend",
            "set_t_low",
            "set_t_high",
            "dust",
            "dust_level",
            "convert",
        ],
    );
    list.extend(excludes(
        "ustat",
        &[
            "checkdup",
            "fa_list",
            "mem",
            "unit",
            "genome_size",
            "sformat",
            "smem",
            "convert",
        ],
    ));
    list.extend(excludes(
        "convert",
        &[
            "checkdup",
            "window",
            "t_extend",
            "t_thres",
            "t_high",
            "t_low",
            "set_t_low",
            "set_t_high",
            "infmt",
            "outfmt",
            "parse_seqids",
            "fa_list",
            "mem",
            "unit",
            "genome_size",
            "dust",
            "dust_level",
            "exclude_ids",
            "ids",
            "text_match",
        ],
    ));
    list
}

// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:79-84
/// windowmasker as NCBI's application runs it (`ncbi_app`): LOSAT's version where NCBI's states its
/// own (`ncbi_version::print`) and no prefix on its diagnostics.
static TOOL: Tool = Tool {
    program: "windowmasker",
    refusal_name: PROGRAM,
    help: HELP,
    descs: DESCS,
    dependencies,
    post_prefix: None,
    components: COMPONENTS,
};

/// The components that `CWinMaskApplication`'s constructor adds: the algorithm, the counts
/// algorithm and the four counts formats, which LOSAT's masker and counts readers implement.
///
/// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:84-95
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker.cpp:57-60,
/// seq_masker_ostat.cpp:41-47 (`STAT_ALGO_*`), seq_masker_ostat_opt_bin.cpp:39-43 (and the other
/// formats' `STAT_FMT_*`)
const COMPONENTS: &[Component] = &[
    component("window-masker-algorithm", "", 0),
    component("windowmasker-statistics-algorithm-version", "", 0),
    component("windowmasker-statistics-format-version", "ascii ", 1),
    component("windowmasker-statistics-format-version", "binary ", 1),
    component("windowmasker-statistics-format-version", "oascii ", 1),
    component("windowmasker-statistics-format-version", "obinary ", 1),
];

const fn component(name: &'static str, prefix: &'static str, minor: u32) -> Component {
    Component {
        name,
        prefix,
        major: 1,
        minor,
        patch: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/app/winmasker/main.cpp
/// Runs `windowmasker` with `words` (the command line after `losat windowmasker`), after its
/// application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, app::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around windowmasker's `Run`:
/// `ncbi_app::run_with`). The library API's `run_windowmasker`.
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, app::run)
}

#[cfg(test)]
mod tests;
