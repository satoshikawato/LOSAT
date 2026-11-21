//! What `CAlnReader::GetSeqEntry` does with the defline of a row (`x_AddMods`, `x_AddTitle`): the
//! `[name=value]` modifiers of the defline (`CTitleParser::Apply`), the handler's choice of what
//! it keeps (`CModHandler::AddMods`), their application to the Bioseq (`CModAdder::Apply` with
//! `CDescrModApply` and `CFeatModApply`), and the title that is left.
//!
//! The Bioseq of NCBI's reader goes to the PSSM of `psiblast -in_msa` (its query), so what
//! matters of the application is:
//!
//! - the messages: warnings are `ERR_POST`s on the standard error stream, an error throws
//!   `CObjReaderParseException` (the reporter has no listener), so the run stops with it;
//! - the title: the text outside the modifiers, then the modifiers that were rejected (the
//!   deprecated ones) and the ones that were skipped (unrecognized, or applied to nothing);
//! - what else the modifiers put on the Bioseq ([`ModEffects`]): descriptors, features and the
//!   Seq-inst. LOSAT does not write them; they show only where the input PSSM is written (the
//!   run refuses there) and the molecule (where the Bioseq stops being protein, the run refuses).

use super::descr::{self, Content};
use crate::ncbi_str::{is_space, string_to_int, string_to_int8, string_to_uint, truncate_spaces};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::TextRun;
use std::collections::BTreeMap;

/// `CModData`: a modifier as the defline gave it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ModData {
    pub(crate) name: Vec<u8>,
    pub(crate) value: Vec<u8>,
}

/// `CSeq_inst::EMol` as the modifiers set it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SetMol {
    Dna,
    Rna,
    Na,
    Aa,
    Other,
}

impl SetMol {
    /// `CBioseq::IsNa()`: `dna`, `rna` and `na`.
    fn is_na(self) -> bool {
        matches!(self, SetMol::Dna | SetMol::Rna | SetMol::Na)
    }
}

/// What the modifiers did to the Bioseq besides its title.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct ModEffects {
    /// The descriptors, the protein feature and the Seq-inst members that the modifiers set.
    pub(crate) content: Content,
    /// The molecule, where a modifier set one (`molecule`, `mol-type`).
    pub(crate) mol: Option<SetMol>,
}

/// The result of the modifiers of one defline.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct ModsApplied {
    /// The title descriptor (`x_AddTitle`): none for a blank title.
    pub(crate) title: Option<Vec<u8>>,
    /// The `ERR_POST(Warning << ..)` messages, in the order NCBI posts them.
    pub(crate) warnings: Vec<Vec<u8>>,
    pub(crate) effects: ModEffects,
}

/// An error that stops the run: the message of the `CObjReaderParseException` and the warnings
/// that NCBI posted before it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ModsFailure {
    pub(crate) warnings: Vec<Vec<u8>>,
    pub(crate) message: Vec<u8>,
}

/// What a modifier of NCBI's reader can do that LOSAT does not follow: a refusal's cause.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ModsStop {
    Failure(ModsFailure),
    /// An exception of NCBI's toolkit that nothing in the reader catches, after the warnings that
    /// NCBI posted before it.
    Exception(Box<crate::ncbi_exception::NcbiException>, Vec<Vec<u8>>),
    /// NCBI's reader does something here that LOSAT does not reproduce.
    Undefined(&'static str),
}

/// `CModHandler::sm_DeprecatedModifiers`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_reader.cpp:92-103
const DEPRECATED: &[&[u8]] = &[
    b"dosage",
    b"transposon-name",
    b"plastid-name",
    b"insertion-seq-name",
    b"old-lineage",
    b"old-name",
    b"gene",
    b"gene-synonym",
    b"allele",
    b"locus-tag",
];

/// `CModHandler::sm_MultipleValuesForbidden`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_reader.cpp:106-128
const SINGLE_VALUE: &[&[u8]] = &[
    b"topology",
    b"molecule",
    b"strand",
    b"gene",
    b"allele",
    b"locus-tag",
    b"protein-desc",
    b"mol-type",
    b"tech",
    b"completeness",
    b"location",
    b"origin",
    b"focus",
    b"taxname",
    b"common",
    b"lineage",
    b"division",
    b"gcode",
    b"mgcode",
    b"pgcode",
];

/// `s_ModNameMap`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_reader.cpp:56-89
const NAME_MAP: &[(&[u8], &[u8])] = &[
    (b"top", b"topology"),
    (b"mol", b"molecule"),
    (b"moltype", b"mol-type"),
    (b"fwd-pcr-primer-name", b"fwd-primer-name"),
    (b"fwd-pcr-primer-names", b"fwd-primer-name"),
    (b"fwd-primer-names", b"fwd-primer-name"),
    (b"fwd-pcr-primer-seq", b"fwd-primer-seq"),
    (b"fwd-pcr-primer-seqs", b"fwd-primer-seq"),
    (b"fwd-primer-seqs", b"fwd-primer-seq"),
    (b"rev-pcr-primer-name", b"rev-primer-name"),
    (b"rev-pcr-primer-names", b"rev-primer-name"),
    (b"rev-primer-names", b"rev-primer-name"),
    (b"rev-pcr-primer-seq", b"rev-primer-seq"),
    (b"rev-pcr-primer-seqs", b"rev-primer-seq"),
    (b"rev-primer-seqs", b"rev-primer-seq"),
    (b"org", b"taxname"),
    (b"organism", b"taxname"),
    (b"div", b"division"),
    (b"notes", b"note"),
    (b"completedness", b"completeness"),
    (b"gene-syn", b"gene-synonym"),
    (b"genesyn", b"gene-synonym"),
    (b"genesynonym", b"gene-synonym"),
    (b"prot", b"protein"),
    (b"prot-desc", b"protein-desc"),
    (b"function", b"activity"),
    (b"secondary", b"secondary-accession"),
    (b"secondary-accessions", b"secondary-accession"),
    (b"keywords", b"keyword"),
    (b"primary", b"primary-accession"),
    (b"primary-accessions", b"primary-accession"),
    (b"projects", b"project"),
    (b"db-xref", b"dbxref"),
    (b"pubmed", b"pmid"),
    (b"ft-url-mod", b"ft-mod"),
    (b"ft-url", b"ft-map"),
    (b"geo-loc-name", b"country"),
];

/// The names of `g_InitModNameSubSrcSubtypeMap`: the `CSubSource` subtypes except the primers
/// (handled elsewhere) and the deprecated ones, with NCBI's extra names.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_to_enum.cpp:121-147
const SUBSOURCE: &[&[u8]] = &[
    b"chromosome",
    b"map",
    b"clone",
    b"subclone",
    b"haplotype",
    b"genotype",
    b"sex",
    b"cell-line",
    b"cell-type",
    b"tissue-type",
    b"clone-lib",
    b"dev-stage",
    b"frequency",
    b"germline",
    b"rearranged",
    b"lab-host",
    b"pop-variant",
    b"tissue-lib",
    b"plasmid-name",
    b"country",
    b"segment",
    b"endogenous-virus-name",
    b"transgenic",
    b"environmental-sample",
    b"isolation-source",
    b"lat-lon",
    b"collection-date",
    b"collected-by",
    b"identified-by",
    b"metagenomic",
    b"mating-type",
    b"linkage-group",
    b"haplogroup",
    b"whole-replicon",
    b"phenotype",
    b"altitude",
    b"other",
    b"sub-clone",
    b"lat-long",
    b"latitude-longitude",
    b"note",
    b"notes",
];

/// The cause of the refusal when a table of this file names a value that the ASN.1 enumeration
/// it is written as does not have (a bug of LOSAT's tables, which the test
/// `the_tables_name_values_of_their_enumerations` keeps out): an error where it could abort.
pub(crate) const INTERNAL_ENUM: &str =
    "a modifier value that LOSAT's table of NCBI's descriptor values cannot write";

/// `name`, when it is a value of the enumeration `member` of the type `ty`.
fn checked(ty: u16, member: &str, name: &'static str) -> Result<&'static str, ModsStop> {
    descr::enum_name(ty, member, name.as_bytes()).ok_or(ModsStop::Undefined(INTERNAL_ENUM))
}

/// The `SubSource` subtype of a name of `SUBSOURCE` (NCBI's extra names map to a type's own).
fn subsource_kind(name: &[u8]) -> Result<&'static str, ModsStop> {
    let kind = match name {
        b"sub-clone" => &b"subclone"[..],
        b"lat-long" | b"latitude-longitude" => b"lat-lon",
        b"note" | b"notes" => b"other",
        other => other,
    };
    descr::enum_name(descr::subsource_type(), "subtype", kind)
        .ok_or(ModsStop::Undefined(INTERNAL_ENUM))
}

/// The `OrgMod` subtype of a name of `ORGMOD`, as `subsource_kind`.
fn orgmod_kind(name: &[u8]) -> Result<&'static str, ModsStop> {
    let kind = match name {
        b"subspecies" => &b"sub-species"[..],
        b"host" | b"specific-host" => b"nat-host",
        other => other,
    };
    descr::enum_name(descr::orgmod_type(), "subtype", kind)
        .ok_or(ModsStop::Undefined(INTERNAL_ENUM))
}

/// `CSubSource::NeedsNoText`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqfeat/SubSource.cpp:231-242
const NEEDS_NO_TEXT: &[&[u8]] = &[
    b"germline",
    b"rearranged",
    b"transgenic",
    b"environmental-sample",
    b"metagenomic",
];

/// The names of `g_InitModNameOrgSubtypeMap`: the `COrgMod` subtypes except the deprecated ones,
/// with NCBI's extra names.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_to_enum.cpp:99-118
const ORGMOD: &[&[u8]] = &[
    b"strain",
    b"substrain",
    b"type",
    b"subtype",
    b"variety",
    b"serotype",
    b"serogroup",
    b"serovar",
    b"cultivar",
    b"pathovar",
    b"chemovar",
    b"biovar",
    b"biotype",
    b"group",
    b"subgroup",
    b"isolate",
    b"common",
    b"acronym",
    b"nat-host",
    b"sub-species",
    b"specimen-voucher",
    b"authority",
    b"forma",
    b"forma-specialis",
    b"ecotype",
    b"synonym",
    b"anamorph",
    b"teleomorph",
    b"breed",
    b"gb-acronym",
    b"gb-anamorph",
    b"gb-synonym",
    b"culture-collection",
    b"bio-material",
    b"metagenome-source",
    b"type-material",
    b"nomenclature",
    b"other",
    b"subspecies",
    b"host",
    b"specific-host",
];

/// `g_BiomolStringToEnum` and what `g_BiomolEnumToMolEnum` makes of each.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_to_enum.cpp:181-232
const BIOMOL: &[(&[u8], &str, SetMol)] = &[
    (b"crna", "cRNA", SetMol::Rna),
    (b"dna", "genomic", SetMol::Dna),
    (b"genomic", "genomic", SetMol::Dna),
    (b"genomicdna", "genomic", SetMol::Dna),
    (b"genomicrna", "genomic", SetMol::Dna),
    (b"mrna", "mRNA", SetMol::Rna),
    (b"ncrna", "ncRNA", SetMol::Rna),
    (b"noncodingrna", "ncRNA", SetMol::Rna),
    (b"othergenetic", "other-genetic", SetMol::Other),
    (b"precursorrna", "pre-RNA", SetMol::Rna),
    (b"ribosomalrna", "rRNA", SetMol::Rna),
    (b"rrna", "rRNA", SetMol::Rna),
    (b"transcribedrna", "transcribed-RNA", SetMol::Rna),
    (b"transfermessengerrna", "tmRNA", SetMol::Rna),
    (b"tmrna", "tmRNA", SetMol::Rna),
    (b"transferrna", "tRNA", SetMol::Rna),
    (b"trna", "tRNA", SetMol::Rna),
];

/// `s_TechStringToEnum`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/descr_mod_apply.cpp:77-104
const TECH: &[(&[u8], &str)] = &[
    (b"?", "unknown"),
    (b"barcode", "barcode"),
    (b"both", "both"),
    (b"compositewgshtgs", "composite-wgs-htgs"),
    (b"concepttrans", "concept-trans"),
    (b"concepttransa", "concept-trans-a"),
    (b"derived", "derived"),
    (b"est", "est"),
    (b"flicdna", "fli-cdna"),
    (b"geneticmap", "genemap"),
    (b"htc", "htc"),
    (b"htgs0", "htgs-0"),
    (b"htgs1", "htgs-1"),
    (b"htgs2", "htgs-2"),
    (b"htgs3", "htgs-3"),
    (b"physicalmap", "physmap"),
    (b"seqpept", "seq-pept"),
    (b"seqpepthomol", "seq-pept-homol"),
    (b"seqpeptoverlap", "seq-pept-overlap"),
    (b"standard", "standard"),
    (b"sts", "sts"),
    (b"survey", "survey"),
    (b"targeted", "targeted"),
    (b"tsa", "tsa"),
    (b"wgs", "wgs"),
];

/// `s_CompletenessStringToEnum`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/descr_mod_apply.cpp:107-116
const COMPLETENESS: &[(&[u8], &str)] = &[
    (b"complete", "complete"),
    (b"hasleft", "has-left"),
    (b"hasright", "has-right"),
    (b"noends", "no-ends"),
    (b"noleft", "no-left"),
    (b"noright", "no-right"),
    (b"partial", "partial"),
];

/// The protein qualifiers of `CFeatModApply::Apply`.
const PROTEIN_QUALS: &[&[u8]] = &[b"protein-desc", b"protein", b"ec-number", b"activity"];

/// `NStr::IsBlank`.
pub(super) fn is_blank(text: &[u8]) -> bool {
    text.iter().all(|&c| is_space(c))
}

/// `NStr::Split(.., delims, fSplit_Tokenize)`: the non-empty runs between delimiters.
pub(super) fn tokenize<'a>(text: &'a [u8], delims: &[u8]) -> Vec<&'a [u8]> {
    text.split(|c| delims.contains(c))
        .filter(|token| !token.is_empty())
        .collect()
}

/// `s_GetNormalizedString` then `s_ModNameMap` (`CModHandler::GetCanonicalName`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_reader.cpp:447-460
fn canonical_name(name: &[u8]) -> Vec<u8> {
    let lower: Vec<u8> = name.iter().map(u8::to_ascii_lowercase).collect();
    let trimmed = truncate_spaces(&lower);
    let is_sep = |c: u8| matches!(c, b'-' | b'_' | b' ');
    let mut normalized: Vec<u8> = Vec::with_capacity(trimmed.len());
    for &c in trimmed {
        if let Some(&last) = normalized.last() {
            if is_sep(last) && is_sep(c) {
                continue;
            }
        }
        normalized.push(c);
    }
    for c in &mut normalized {
        if matches!(*c, b'_' | b' ') {
            *c = b'-';
        }
    }
    match NAME_MAP
        .iter()
        .find(|(from, _)| *from == normalized.as_slice())
    {
        Some((_, to)) => to.to_vec(),
        None => normalized,
    }
}

/// `g_GetNormalizedModVal`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_to_enum.cpp:35-64
pub(super) fn normalized_value(value: &[u8]) -> Vec<u8> {
    let edge = |c: &u8| matches!(c, b' ' | b'\t' | b'"' | b'\'' | b'-' | b'_');
    let start = value.iter().position(|c| !edge(c));
    let trimmed: &[u8] = match start {
        Some(start) => {
            let end = value.iter().rposition(|c| !edge(c)).unwrap_or(start);
            &value[start..=end]
        }
        None => value,
    };
    trimmed
        .iter()
        .filter(|&&c| !matches!(c, b'-' | b'_' | b':') && !is_space(c))
        .map(u8::to_ascii_lowercase)
        .collect()
}

/// `CTitleParser::Apply`: the modifiers of the title and the text around them.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_reader.cpp:704-748
fn title_parser_apply(title: &[u8]) -> (Vec<ModData>, Vec<u8>) {
    let mut mods = Vec::new();
    let mut remainder: Vec<u8> = Vec::new();
    let append = |remainder: &mut Vec<u8>, text: &[u8]| {
        let text = truncate_spaces(text);
        if !text.is_empty() {
            if !remainder.is_empty() {
                remainder.push(b' ');
            }
            remainder.extend_from_slice(text);
        }
    };
    let mut start_pos = 0;
    while start_pos < title.len() {
        let mut lb_pos = start_pos;
        match crate::blastinput::title_parser::find_brackets(title, &mut lb_pos) {
            Some((eq_pos, end_pos)) => {
                if lb_pos > start_pos {
                    append(&mut remainder, &title[start_pos..lb_pos]);
                }
                mods.push(ModData {
                    name: truncate_spaces(&title[lb_pos + 1..eq_pos]).to_vec(),
                    value: truncate_spaces(&title[eq_pos + 1..end_pos]).to_vec(),
                });
                start_pos = end_pos + 1;
            }
            None => {
                append(&mut remainder, &title[start_pos..]);
                return (mods, remainder);
            }
        }
    }
    (mods, remainder)
}

type ModMap = BTreeMap<Vec<u8>, Vec<ModData>>;

fn message(parts: &[&[u8]]) -> Vec<u8> {
    parts.concat()
}

/// `CModHandler::AddMods` (`eAppendReplace`, no excluded or ignored modifiers, the reporter of
/// `CDefaultModErrorReporter`): the modifiers it keeps by canonical name, and the ones it
/// rejects (deprecated). Warnings go to `warnings`; an error stops.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_reader.cpp:231-326
fn add_mods(
    mods: &[ModData],
    warnings: &mut Vec<Vec<u8>>,
) -> Result<(ModMap, Vec<ModData>), Vec<u8>> {
    let mut rejected = Vec::new();
    let mut seen: Vec<Vec<u8>> = Vec::new();
    let mut accepted: ModMap = BTreeMap::new();
    for modifier in mods {
        let canonical = canonical_name(&modifier.name);
        let multiple_allowed = !SINGLE_VALUE.contains(&canonical.as_slice());
        if DEPRECATED.contains(&canonical.as_slice()) {
            warnings.push(message(&[
                b"Use of the following modifier in a sequence file is discouraged and the information will be ignored: ",
                &modifier.name,
                b".",
            ]));
            rejected.push(modifier.clone());
            continue;
        }
        if seen.contains(&canonical) {
            let same_value = accepted
                .get(&canonical)
                .and_then(|list| list.first())
                .is_some_and(|first| first.value.eq_ignore_ascii_case(&modifier.value));
            if same_value {
                warnings.push(message(&[
                    b"Duplicated modifier value detected, ignoring duplicate, no action required: ",
                    &modifier.name,
                    b"=",
                    &modifier.value,
                    b".",
                ]));
                continue;
            }
            if !multiple_allowed {
                return Err(message(&[
                    b"Conflicting modifiers detected. Provide one modifier with one value for: ",
                    &modifier.name,
                    b".",
                ]));
            }
        } else {
            seen.push(canonical.clone());
        }
        accepted
            .entry(canonical)
            .or_default()
            .push(modifier.clone());
    }
    Ok((accepted, rejected))
}

/// The state of `CModAdder::Apply` over one Bioseq.
struct Adder {
    /// `CSeq_inst::mol` as the modifiers have set it so far.
    mol: SetMol,
    mol_set: bool,
    effects: ModEffects,
    warnings: Vec<Vec<u8>>,
    skipped: Vec<ModData>,
}

fn invalid_value(modifier: &ModData, extra: &[u8]) -> Vec<u8> {
    let mut text = message(&[
        b"Invalid value: ",
        &modifier.name,
        b"=",
        &modifier.value,
        b".",
    ]);
    if !extra.is_empty() {
        text.push(b' ');
        text.extend_from_slice(extra);
    }
    text
}

impl Adder {
    fn single_value(values: &[ModData]) -> &[u8] {
        &values[0].value
    }

    /// `CDescrModApply::Apply`: whether the modifier is one of its own (`Ok(true)`), where an
    /// invalid value is an error.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/descr_mod_apply.cpp:186-232
    fn descr_apply(&mut self, canonical: &[u8], values: &[ModData]) -> Result<bool, ModsStop> {
        let fail = |modifier: &ModData, extra: &[u8]| {
            Err(ModsStop::Failure(ModsFailure {
                warnings: Vec::new(),
                message: invalid_value(modifier, extra),
            }))
        };
        let value = Self::single_value(values);
        let texts = || -> Vec<Vec<u8>> { values.iter().map(|m| m.value.clone()).collect() };
        let split_all = |delims: &[u8]| -> Vec<Vec<u8>> {
            values
                .iter()
                .flat_map(|modifier| tokenize(&modifier.value, delims))
                .map(<[u8]>::to_vec)
                .collect()
        };
        let content = &mut self.effects.content;
        match canonical {
            b"location" => {
                let key = normalized_value(value);
                let Some(name) = descr::biosource_enum(
                    "genome",
                    &key,
                    &[
                        (b"mitochondrial", "mitochondrion"),
                        (b"provirus", "proviral"),
                        (b"extrachromosomal", "extrachrom"),
                        (b"insertion sequence", "insertion-seq"),
                    ],
                ) else {
                    return fail(&values[0], b"");
                };
                content.set_genome(checked(named::BIOSOURCE, "genome", name)?);
            }
            b"origin" => {
                let key = normalized_value(value);
                let Some(name) = descr::biosource_enum(
                    "origin",
                    &key,
                    &[(b"natural mutant", "natmut"), (b"mutant", "mut")],
                ) else {
                    return fail(&values[0], b"");
                };
                content.set_origin(checked(named::BIOSOURCE, "origin", name)?);
            }
            b"focus" => {
                if value.eq_ignore_ascii_case(b"true") {
                    content.set_focus();
                } else if value.eq_ignore_ascii_case(b"false") {
                    return fail(&values[0], b"");
                }
            }
            b"fwd-primer-name" => content.set_primers(true, true, &texts()),
            b"fwd-primer-seq" => content.set_primers(true, false, &texts()),
            b"rev-primer-name" => content.set_primers(false, true, &texts()),
            b"rev-primer-seq" => content.set_primers(false, false, &texts()),
            b"taxname" => content.set_taxname(value),
            b"common" => content.set_common(value),
            b"dbxref" => content.set_dbxref(&texts()),
            b"lineage" => content.set_lineage(value),
            b"division" => content.set_division(value),
            b"taxid" => {
                let Ok(id) = string_to_int(value, 0) else {
                    return fail(&values[0], b"Integer value expected.");
                };
                content.set_taxid(i64::from(id));
            }
            b"gcode" | b"mgcode" | b"pgcode" => {
                let Ok(code) = string_to_int(value, 0) else {
                    return fail(&values[0], b"Integer value expected.");
                };
                content.set_code(canonical, i64::from(code));
            }
            b"sra" | b"bioproject" | b"biosample" => {
                let label: &[u8] = match canonical {
                    b"sra" => b"Sequence Read Archive",
                    b"biosample" => b"BioSample",
                    _ => b"BioProject",
                };
                let all = split_all(b",; \t");
                if !all.is_empty() {
                    content.set_dblink(label, all);
                }
            }
            b"keyword" => {
                let all = split_all(b",; \t");
                if !all.is_empty() {
                    content.set_keywords(all);
                }
            }
            b"mol-type" => {
                let key = normalized_value(value);
                let Some((_, biomol, mol)) = BIOMOL.iter().find(|(k, _, _)| *k == key.as_slice())
                else {
                    return fail(&values[0], b"");
                };
                content.set_biomol(checked(named::MOLINFO, "biomol", biomol)?);
                self.mol = *mol;
                self.mol_set = true;
            }
            b"completeness" => {
                let key = normalized_value(value);
                let Some((_, name)) = COMPLETENESS.iter().find(|(k, _)| *k == key.as_slice())
                else {
                    return fail(&values[0], b"");
                };
                content.set_completeness(checked(named::MOLINFO, "completeness", name)?);
            }
            b"tech" => {
                let key = normalized_value(value);
                let Some((_, name)) = TECH.iter().find(|(k, _)| *k == key.as_slice()) else {
                    return fail(&values[0], b"");
                };
                content.set_tech(checked(named::MOLINFO, "tech", name)?);
            }
            b"primary-accession" => {
                let mut accessions: Vec<Vec<u8>> = Vec::new();
                for modifier in values {
                    // `CUtf8::AsUTF8` of a token that is not UTF-8 throws, and the modifier is
                    // reported as an invalid value.
                    let tokens = tokenize(&modifier.value, b",; \t");
                    if tokens
                        .iter()
                        .any(|token| std::str::from_utf8(token).is_err())
                    {
                        return fail(modifier, b"");
                    }
                    accessions.extend(tokens.into_iter().map(<[u8]>::to_vec));
                }
                if !accessions.is_empty() {
                    content.set_tpa(accessions);
                }
            }
            b"secondary-accession" => {
                let ids = secondary_ids(values)?;
                // `CModAdder::x_SetHist`: a Seq-id of each (`CSeq_id(string)`), whose failure
                // is an exception that nothing in the reader catches. The ids are made one at a
                // time here to be checked, and again when the descriptor is written.
                let mut any = false;
                for run in &ids {
                    let mut failure = None;
                    run.for_each(&mut |id| {
                        any = true;
                        if failure.is_none() {
                            if let Err(error) = crate::seqid::parse::seq_id_from_text(
                                id,
                                crate::seqid::accguide::parse_flags::ANY_RAW,
                                &mut |_| {},
                            ) {
                                failure = Some(error);
                            }
                        }
                    });
                    if let Some(error) = failure {
                        return Err(ModsStop::Exception(Box::new(error.exception()), Vec::new()));
                    }
                }
                if any {
                    content.replaces = Some(ids.clone());
                }
                content.set_extra_accessions(ids);
            }
            b"project" => {
                let mut ids: Vec<i64> = Vec::new();
                for modifier in values {
                    for token in tokenize(&modifier.value, b",; \t") {
                        match string_to_uint(token) {
                            Ok(id) => ids.push(i64::from(id)),
                            Err(_) => return fail(modifier, b""),
                        }
                    }
                }
                if !ids.is_empty() {
                    content.set_genome_projects(ids);
                }
            }
            b"comment" => {
                for modifier in values {
                    content.add_comment(&modifier.value);
                }
            }
            b"ft-map" | b"ft-mod" => {
                let label: &[u8] = if canonical == b"ft-map" {
                    b"Map-FileTrackURL"
                } else {
                    b"BaseModification-FileTrackURL"
                };
                for modifier in values {
                    content.add_file_track(label, &modifier.value);
                }
            }
            b"pmid" => {
                for modifier in values {
                    match string_to_int8(&modifier.value, 0) {
                        Ok(pmid) => content.add_pmid(pmid),
                        Err(_) => return fail(&values[0], b"Expected integer value."),
                    }
                }
            }
            name if SUBSOURCE.contains(&name) => {
                let kind = subsource_kind(name)?;
                if NEEDS_NO_TEXT.contains(&kind.as_bytes()) {
                    for modifier in values {
                        if !modifier.value.eq_ignore_ascii_case(b"true") {
                            return fail(modifier, b"");
                        }
                    }
                }
                content.add_subsources(kind, &texts());
            }
            name if ORGMOD.contains(&name) => {
                let kind = orgmod_kind(name)?;
                content.add_orgmods(kind, &texts());
            }
            _ => return Ok(false),
        }
        Ok(true)
    }

    /// `CModAdder::x_TrySeqInstMod`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_reader.cpp:549-575
    fn seq_inst_apply(&mut self, canonical: &[u8], values: &[ModData]) -> Result<bool, ModsStop> {
        let normalized = normalized_value(Self::single_value(values));
        let fail = || {
            Err(ModsStop::Failure(ModsFailure {
                warnings: Vec::new(),
                message: invalid_value(&values[0], b""),
            }))
        };
        match canonical {
            b"strand" => {
                self.effects.content.strand = Some(match normalized.as_slice() {
                    b"single" => "ss",
                    b"double" => "ds",
                    b"mixed" => "mixed",
                    b"other" => "other",
                    _ => return fail(),
                });
            }
            b"molecule" => {
                let mol = match normalized.as_slice() {
                    b"dna" => SetMol::Dna,
                    b"rna" => SetMol::Rna,
                    b"aa" => SetMol::Aa,
                    b"na" => SetMol::Na,
                    b"other" => SetMol::Other,
                    _ => return fail(),
                };
                self.mol = mol;
                self.mol_set = true;
            }
            b"topology" => {
                self.effects.content.topology = Some(match normalized.as_slice() {
                    b"linear" => "linear",
                    b"circular" => "circular",
                    b"tandem" => "tandem",
                    b"other" => "other",
                    _ => return fail(),
                });
            }
            _ => return Ok(false),
        }
        Ok(true)
    }

    /// `CFeatModApply::Apply`: the protein qualifiers (a warning on a nucleotide Bioseq).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/feature_mod_apply.cpp:62-106
    fn feat_apply(&mut self, canonical: &[u8], values: &[ModData]) -> bool {
        if !PROTEIN_QUALS.contains(&canonical) {
            return false;
        }
        if self.mol.is_na() {
            self.warnings.push(message(&[
                b"Cannot apply protein modifier to nucleotide sequence. The following modifier will be ignored: ",
                &values[0].name,
            ]));
            self.skipped.extend(values.iter().cloned());
        } else {
            let prot = self.effects.content.protein();
            let texts = || -> Vec<Vec<u8>> { values.iter().map(|m| m.value.clone()).collect() };
            match canonical {
                b"protein-desc" => prot.desc = Some(Self::single_value(values).to_vec()),
                b"protein" => prot.name = Some(texts()),
                b"ec-number" => prot.ec = Some(texts()),
                _ => prot.activity = Some(texts()),
            }
        }
        true
    }
}

/// The cause of the refusal of a secondary accession whose range runs backwards by more than
/// one: `std::list::insert` over the range's iterators then counts up from the range's start
/// past `INT_MAX` (signed overflow), and NCBI BLAST+ 2.17.0 runs until its memory or time is
/// gone (probe `f6b.m2.028`: no answer in 280 s).
pub(crate) const BACKWARD_RANGE: &str =
    "NCBI BLAST+ 2.17.0 reads a secondary-accession range that runs backwards (AB000003-AB000001) as a loop past the largest integer, which does not end";

/// The cause of the refusal of a secondary-accession range whose numbers do not fit 64 bits.
pub(crate) const RANGE_OVERFLOW: &str =
    "a secondary-accession range with a number of more than 18 digits, which NCBI BLAST+ 2.17.0 reads into a 32-bit integer (undefined behavior) and LOSAT cannot expand";

/// `SSeqIdRange` of one value as a run of its accessions: `None` where the constructor throws (the
/// value is then taken as it is). NCBI reads the numbers into a signed 32-bit `int`, which wraps
/// for ten digits or more (signed overflow: undefined behavior, Owner 2026-10-22 「広い型になおし
/// ましょう」, PD-LOSAT-NCBI-DEFECTS): LOSAT reads them into 64 bits and gives the accessions with
/// the digits as written, as NCBI does for every number that fits.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:3602-3642
fn seq_id_range(value: &[u8]) -> Result<Option<TextRun>, ModsStop> {
    let alpha = |c: u8| c.is_ascii_alphabetic();
    let n = value.len();
    let mut pos = 0;
    while pos < n && alpha(value[pos]) {
        pos += 1;
    }
    let prefix = &value[..pos];
    let mut start: i64 = 0;
    let mut digits = 0usize;
    let append = |number: i64, digit: u8| -> Result<i64, ModsStop> {
        number
            .checked_mul(10)
            .and_then(|number| number.checked_add(i64::from(digit - b'0')))
            .ok_or(ModsStop::Undefined(RANGE_OVERFLOW))
    };
    while pos < n && value[pos].is_ascii_digit() {
        start = append(start, value[pos])?;
        digits += 1;
        pos += 1;
    }
    let stop;
    if pos == n {
        stop = start;
    } else {
        if value[pos] != b'-' {
            return Ok(None);
        }
        pos += 1;
        let second = pos;
        while pos < n && alpha(value[pos]) {
            pos += 1;
        }
        let second = &value[second..pos];
        if !second.is_empty() && second != prefix {
            return Ok(None);
        }
        if pos + digits != n {
            return Ok(None);
        }
        let mut number: i64 = 0;
        while pos < n && value[pos].is_ascii_digit() {
            number = append(number, value[pos])?;
            pos += 1;
        }
        stop = number;
    }
    if start - 1 > stop {
        return Err(ModsStop::Undefined(BACKWARD_RANGE));
    }
    Ok(Some(TextRun::Numbered {
        prefix: prefix.to_vec(),
        first: start,
        count: (stop - start).saturating_add(1),
        digits,
    }))
}

/// The ids of a `secondary-accession` modifier (`x_SetGBblockIds`, `x_SetHist`): the values split
/// at `,; \t`, each trimmed, a range as one run (written when the descriptor is, not held).
fn secondary_ids(values: &[ModData]) -> Result<Vec<TextRun>, ModsStop> {
    let mut ids = Vec::new();
    for modifier in values {
        for token in tokenize(&modifier.value, b",; \t") {
            let value = truncate_spaces(token);
            match seq_id_range(value)? {
                Some(range) => ids.push(range),
                None => ids.push(TextRun::One(value.to_vec())),
            }
        }
    }
    Ok(ids)
}

/// `CAlnReader::x_AddMods` of one defline: its title and what the modifiers did.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_reader.cpp:807-849
pub(crate) fn apply_defline(defline: &[u8], protein: bool) -> Result<ModsApplied, ModsStop> {
    let mut applied = ModsApplied::default();
    if is_blank(defline) {
        return Ok(applied);
    }
    let (mods, mut remainder) = title_parser_apply(defline);
    if mods.is_empty() && is_blank(&remainder) {
        return Ok(applied);
    }
    let mut warnings = Vec::new();
    let (accepted, rejected) = match add_mods(&mods, &mut warnings) {
        Ok(result) => result,
        Err(message) => {
            return Err(ModsStop::Failure(ModsFailure { warnings, message }));
        }
    };
    let mut adder = Adder {
        mol: if protein { SetMol::Aa } else { SetMol::Dna },
        mol_set: false,
        effects: ModEffects::default(),
        warnings,
        skipped: Vec::new(),
    };
    for (canonical, values) in &accepted {
        let step = (|| -> Result<(), ModsStop> {
            if adder.descr_apply(canonical, values)? {
                return Ok(());
            }
            if adder.seq_inst_apply(canonical, values)? || adder.feat_apply(canonical, values) {
                return Ok(());
            }
            for modifier in values {
                adder.skipped.push(modifier.clone());
                adder
                    .warnings
                    .push(message(&[b"Unrecognized modifier: ", &modifier.name, b"."]));
            }
            Ok(())
        })();
        if let Err(stop) = step {
            return Err(match stop {
                ModsStop::Failure(failure) => ModsStop::Failure(ModsFailure {
                    warnings: adder.warnings,
                    message: failure.message,
                }),
                ModsStop::Exception(exception, _) => ModsStop::Exception(exception, adder.warnings),
                other => other,
            });
        }
    }
    for modifier in rejected.iter().chain(adder.skipped.iter()) {
        remainder.extend_from_slice(b" [");
        remainder.extend_from_slice(&modifier.name);
        remainder.push(b'=');
        remainder.extend_from_slice(&modifier.value);
        remainder.push(b']');
    }
    let title = truncate_spaces(&remainder);
    if !is_blank(title) {
        applied.title = Some(title.to_vec());
    }
    applied.warnings = adder.warnings;
    applied.effects = adder.effects;
    if adder.mol_set {
        applied.effects.mol = Some(adder.mol);
    }
    Ok(applied)
}

#[cfg(test)]
#[path = "mods_tests.rs"]
mod tests;
