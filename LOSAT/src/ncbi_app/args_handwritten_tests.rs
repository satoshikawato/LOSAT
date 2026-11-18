//! Temporary (line RG DECISIONS K4i-1, test (b)): the descriptions that NCBI's XML gives
//! (`Arguments`) are the hand-written ones of dev f8dd8bed2 (copied below verbatim, without
//! their comments): the same keys and flags, the hand-written order is that of the usage
//! (`usage_position`), and the `SetDependency` calls, made into NCBI's `multimap`, are the XML's
//! dependencies. Removed with the hand-written lists once it passes.

use super::*;

const fn key(
    name: &'static str,
    kind: Kind,
    default: Option<&'static str>,
    constraint: Option<Constraint>,
) -> Desc {
    Desc {
        name,
        kind,
        default,
        constraint,
        mandatory: false,
    }
}

const fn mandatory_key(name: &'static str, kind: Kind) -> Desc {
    Desc {
        name,
        kind,
        default: None,
        constraint: None,
        mandatory: true,
    }
}

mod blastdbcmd {
    use super::*;

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
        key("ipg", Kind::Integer, None, Some(Constraint::AtLeast(0.0))),
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
            Some(Constraint::AtLeast(1.0)),
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
}

mod dustmasker {
    use super::*;

    pub(super) const DESCS: &[Desc] = &[
        key("in", Kind::InputFile, Some("-"), None),
        key("out", Kind::OutputFile, Some("-"), None),
        key("window", Kind::Integer, Some("64"), None),
        key("level", Kind::Integer, Some("20"), None),
        key("linker", Kind::Integer, Some("1"), None),
        key("infmt", Kind::Text, Some("fasta"), None),
        key(
            "outfmt",
            Kind::Text,
            Some("interval"),
            Some(Constraint::Strings(&[
                "acclist",
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
        flag("parse_seqids"),
        flag("hard_masking"),
    ];

    pub(super) fn dependencies() -> Vec<(&'static str, Dep, &'static str)> {
        Vec::new()
    }
}

mod segmasker {
    use super::*;

    pub(super) const DESCS: &[Desc] = &[
        key("in", Kind::InputFile, Some("-"), None),
        key("out", Kind::OutputFile, Some("-"), None),
        key(
            "infmt",
            Kind::Text,
            Some("fasta"),
            Some(Constraint::Strings(&["blastdb", "fasta"])),
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
        key("window", Kind::Integer, Some("12"), None),
        key("locut", Kind::Double, Some("2.2"), None),
        key("hicut", Kind::Double, Some("2.5"), None),
    ];

    pub(super) fn dependencies() -> Vec<(&'static str, Dep, &'static str)> {
        Vec::new()
    }
}

mod convert2blastmask {
    use super::*;

    pub(super) const DESCS: &[Desc] = &[
        key("in", Kind::InputFile, Some("-"), None),
        key("out", Kind::OutputFile, Some("-"), None),
        key(
            "outfmt",
            Kind::Text,
            Some("maskinfo_asn1_text"),
            Some(Constraint::Strings(&[
                "interval",
                "maskinfo_asn1_bin",
                "maskinfo_asn1_text",
                "maskinfo_xml",
            ])),
        ),
        flag("parse_seqids"),
        mandatory_key("masking_algorithm", Kind::Text),
        mandatory_key("masking_options", Kind::Text),
    ];

    pub(super) fn dependencies() -> Vec<(&'static str, Dep, &'static str)> {
        Vec::new()
    }
}

/// The hand-written `SetDependency` calls as NCBI's `multimap` holds them (dev's
/// `args::dependencies`), its `eRequires` entries first, as the XML prints them.
fn multimap(
    calls: Vec<(&'static str, Dep, &'static str)>,
) -> Vec<(&'static str, Dep, &'static str)> {
    let mut deps = Vec::new();
    for (first, dep, second) in calls {
        deps.push((first, dep, second));
        if dep == Dep::Excludes {
            deps.push((second, dep, first));
        }
    }
    deps.sort_by_key(|&(first, _, _)| first);
    let (mut requires, excludes): (Vec<_>, Vec<_>) = deps
        .into_iter()
        .partition(|&(_, dep, _)| dep == Dep::Requires);
    requires.extend(excludes);
    requires
}

#[test]
fn the_xml_gives_the_hand_written_descriptions() {
    type Calls = fn() -> Vec<(&'static str, Dep, &'static str)>;
    let programs: [(&str, &str, &[Desc], Calls); 4] = [
        (
            "blastdbcmd",
            include_str!("../blastdbcmd/help.txt"),
            blastdbcmd::DESCS,
            blastdbcmd::dependencies,
        ),
        (
            "dustmasker",
            include_str!("../dustmasker/help.txt"),
            dustmasker::DESCS,
            dustmasker::dependencies,
        ),
        (
            "segmasker",
            include_str!("../segmasker/help.txt"),
            segmasker::DESCS,
            segmasker::dependencies,
        ),
        (
            "convert2blastmask",
            include_str!("../convert2blastmask/help.txt"),
            convert2blastmask::DESCS,
            convert2blastmask::dependencies,
        ),
    ];
    for (program, help, hand_written, calls) in programs {
        let arguments = arguments(program).unwrap();
        assert_eq!(arguments.descs.len(), hand_written.len(), "{program}");
        for desc in hand_written {
            let generated = arguments.descs.iter().find(|g| g.name == desc.name);
            assert_eq!(generated, Some(desc), "{program}");
        }
        let positions: Vec<usize> = hand_written
            .iter()
            .map(|desc| usage_position(help, desc.name))
            .collect();
        assert!(
            positions.windows(2).all(|pair| pair[0] < pair[1]) && !positions.contains(&usize::MAX),
            "{program}: {positions:?}"
        );
        assert!(arguments.aliases.is_empty(), "{program}");
        assert_eq!(arguments.dependencies, &multimap(calls())[..], "{program}");
    }
}
