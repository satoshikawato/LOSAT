//! The arguments of `losat blastdbcmd`: the descriptions of blastdbcmd 2.17.0
//! (`CBlastDBCmdApp::Init`) that NCBI's argument reader reads (`ncbi_app::args`), so that every
//! command line that NCBI's reader refuses is refused with NCBI's text, also for the options that
//! LOSAT's BLASTDBCMD does not serve (the program refuses those after the reader accepts them).
//! With them `string_to_int8_message`, `NStr::StringToInt8` (`ncbi_str`) for the values
//! blastdbcmd converts itself (`-mask_sequence_with`, the lines of its batch files).

use crate::ncbi_app::args::{flag, key, Constraint, Dep, Desc, Kind};

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:901-919
/// The arguments of `CBlastDBCmdApp::Init` (blastdbcmd.cpp:901-1167) in the order it describes
/// them; NCBI's reader adds the help flags before them and `-logfile` and `-conffile` after them
/// (`ncbi_app::args`). The file keys are `eInputFile` (`-entry_batch`, `-ipg_batch`,
/// `-taxidlist`) and `eOutputFile` (`-out`; blastdbcmd.cpp:928-931, 945-948, 1029-1030,
/// 1042-1043).
pub(super) const DESCS: &[Desc] = &[
    key("db", Kind::Text, Some("nr"), None),
    key(
        "dbtype",
        Kind::Text,
        Some("guess"),
        Some(Constraint::Strings(&["guess", "nucl", "prot"])),
    ),
    key("entry", Kind::Text, None, None),
    key("entry_batch", Kind::InputFile, None, None),
    key("ipg", Kind::Integer, None, Some(Constraint::AtLeast(0))),
    key("ipg_batch", Kind::InputFile, None, None),
    flag("info"),
    flag("metadata"),
    key("metadata_output_prefix", Kind::Text, None, None),
    flag("tax_info"),
    key("range", Kind::Text, None, None),
    key(
        "strand",
        Kind::Text,
        Some("plus"),
        Some(Constraint::Strings(&["minus", "plus"])),
    ),
    key("mask_sequence_with", Kind::Text, None, None),
    key("taxids", Kind::Text, None, None),
    key("taxidlist", Kind::InputFile, None, None),
    flag("no_taxid_expansion"),
    key("out", Kind::OutputFile, Some("-"), None),
    key("outfmt", Kind::Text, Some("%f"), None),
    flag("target_only"),
    flag("get_dups"),
    key(
        "line_length",
        Kind::Integer,
        Some("80"),
        Some(Constraint::AtLeast(1)),
    ),
    flag("ctrl_a"),
    flag("show_blastdb_search_path"),
    key("list", Kind::Text, None, None),
    flag("remove_redundant_dbs"),
    flag("recursive"),
    key("list_outfmt", Kind::Text, Some("%f %p"), None),
    flag("exact_length"),
    flag("long_seqids"),
];

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:955-961
// `pig` and `mask_sequence` name no argument, so they never apply.
const INFO_EXCLUSIONS: &[&str] = &[
    "entry",
    "entry_batch",
    "outfmt",
    "strand",
    "target_only",
    "ctrl_a",
    "get_dups",
    "pig",
    "range",
    "mask_sequence",
    "list",
    "remove_redundant_dbs",
    "recursive",
    "list_outfmt",
    "metadata",
    "metadata_output_prefix",
    "taxidlist",
    "taxids",
    "no_taxid_expansion",
];

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:969-977
const METADATA_EXCLUSIONS: &[&str] = &[
    "entry",
    "entry_batch",
    "outfmt",
    "strand",
    "target_only",
    "ctrl_a",
    "get_dups",
    "pig",
    "range",
    "mask_sequence",
    "list",
    "remove_redundant_dbs",
    "recursive",
    "list_outfmt",
    "info",
    "taxidlist",
    "taxids",
    "no_taxid_expansion",
];

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:994-998
const TAX_INFO_EXCLUSIONS: &[&str] = &[
    "info",
    "entry",
    "entry_batch",
    "strand",
    "target_only",
    "ctrl_a",
    "get_dups",
    "pig",
    "range",
    "mask_sequence",
    "list",
    "remove_redundant_dbs",
    "recursive",
    "list_outfmt",
    "taxidlist",
    "taxids",
    "no_taxid_expansion",
];

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1111-1113
const DISCOVERY_EXCLUSIONS: &[&str] = &[
    "entry",
    "entry_batch",
    "outfmt",
    "strand",
    "target_only",
    "ctrl_a",
    "get_dups",
    "pig",
    "range",
    "db",
    "info",
    "mask_sequence",
    "line_length",
];

/// The `SetDependency` calls of blastdbcmd (blastdbcmd.cpp:901-1167), in their order; NCBI's
/// reader holds them in its `multimap` (`ncbi_app::args::dependencies`).
pub(super) fn dependencies() -> Vec<(&'static str, Dep, &'static str)> {
    let mut deps = Vec::new();
    let mut set = |first: &'static str, dep: Dep, second: &'static str| {
        deps.push((first, dep, second));
    };
    for second in ["entry", "range", "strand", "mask_sequence_with"] {
        set("entry_batch", Dep::Excludes, second);
    }
    for second in ["entry", "entry_batch", "target_only", "ipg_batch"] {
        set("ipg", Dep::Excludes, second);
    }
    for second in [
        "entry",
        "entry_batch",
        "range",
        "strand",
        "mask_sequence_with",
    ] {
        set("ipg_batch", Dep::Excludes, second);
    }
    for &second in INFO_EXCLUSIONS {
        set("info", Dep::Excludes, second);
    }
    for &second in METADATA_EXCLUSIONS {
        set("metadata", Dep::Excludes, second);
    }
    set("metadata_output_prefix", Dep::Requires, "metadata");
    for &second in TAX_INFO_EXCLUSIONS {
        set("tax_info", Dep::Excludes, second);
    }
    for second in ["entry", "entry_batch", "pig"] {
        set("taxids", Dep::Excludes, second);
    }
    for second in ["entry", "entry_batch", "pig", "taxids"] {
        set("taxidlist", Dep::Excludes, second);
    }
    for second in ["entry", "entry_batch", "pig"] {
        set("no_taxid_expansion", Dep::Excludes, second);
    }
    set("get_dups", Dep::Excludes, "target_only");
    for &second in DISCOVERY_EXCLUSIONS {
        for first in [
            "list",
            "recursive",
            "remove_redundant_dbs",
            "list_outfmt",
            "show_blastdb_search_path",
        ] {
            set(first, Dep::Excludes, second);
        }
    }
    for second in ["list", "recursive", "list_outfmt", "remove_redundant_dbs"] {
        set("show_blastdb_search_path", Dep::Excludes, second);
    }
    set("exact_length", Dep::Requires, "info");
    set("long_seqids", Dep::Excludes, "info");
    deps
}

/// `NStr::StringToInt8(value)` (`ncbi_str::string_to_int8`): the number, or the message of
/// NCBI's exception as `GetMsg` gives it (without the position).
pub(super) fn string_to_int8_message(value: &str) -> Result<i64, String> {
    crate::ncbi_str::string_to_int8(value.as_bytes(), 0)
        .map_err(|error| String::from_utf8_lossy(error.get_msg()).into_owned())
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod tests;
