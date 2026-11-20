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

use crate::ncbi_str::{string_to_int, string_to_int8, string_to_uint, truncate_spaces};
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
    /// A descriptor (BioSource, MolInfo, GB-block, user object, comment, publication), a feature
    /// (protein) or a Seq-inst member (strand, topology, history) that the modifiers set.
    pub(crate) non_title: bool,
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

/// `BioSource.genome` values by `g_GetNormalizedModVal`.
const GENOME: &[&[u8]] = &[
    b"unknown",
    b"genomic",
    b"chloroplast",
    b"chromoplast",
    b"kinetoplast",
    b"mitochondrion",
    b"plastid",
    b"macronuclear",
    b"extrachrom",
    b"plasmid",
    b"transposon",
    b"insertionseq",
    b"cyanelle",
    b"proviral",
    b"virion",
    b"nucleomorph",
    b"apicoplast",
    b"leucoplast",
    b"proplastid",
    b"endogenousvirus",
    b"hydrogenosome",
    b"chromosome",
    b"chromatophore",
    b"plasmidinmitochondrion",
    b"plasmidinplastid",
    b"mitochondrial",
    b"provirus",
    b"extrachromosomal",
    b"insertionsequence",
];

/// `BioSource.origin` values by `g_GetNormalizedModVal`.
const ORIGIN: &[&[u8]] = &[
    b"unknown",
    b"natural",
    b"natmut",
    b"mut",
    b"artificial",
    b"synthetic",
    b"other",
    b"naturalmutant",
    b"mutant",
];

/// `g_BiomolStringToEnum` and what `g_BiomolEnumToMolEnum` makes of each.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_to_enum.cpp:181-232
const BIOMOL: &[(&[u8], SetMol)] = &[
    (b"crna", SetMol::Rna),
    (b"dna", SetMol::Dna),
    (b"genomic", SetMol::Dna),
    (b"genomicdna", SetMol::Dna),
    (b"genomicrna", SetMol::Dna),
    (b"mrna", SetMol::Rna),
    (b"ncrna", SetMol::Rna),
    (b"noncodingrna", SetMol::Rna),
    (b"othergenetic", SetMol::Other),
    (b"precursorrna", SetMol::Rna),
    (b"ribosomalrna", SetMol::Rna),
    (b"rrna", SetMol::Rna),
    (b"transcribedrna", SetMol::Rna),
    (b"transfermessengerrna", SetMol::Rna),
    (b"tmrna", SetMol::Rna),
    (b"transferrna", SetMol::Rna),
    (b"trna", SetMol::Rna),
];

/// `s_TechStringToEnum`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/descr_mod_apply.cpp:77-104
const TECH: &[&[u8]] = &[
    b"?",
    b"barcode",
    b"both",
    b"compositewgshtgs",
    b"concepttrans",
    b"concepttransa",
    b"derived",
    b"est",
    b"flicdna",
    b"geneticmap",
    b"htc",
    b"htgs0",
    b"htgs1",
    b"htgs2",
    b"htgs3",
    b"physicalmap",
    b"seqpept",
    b"seqpepthomol",
    b"seqpeptoverlap",
    b"standard",
    b"sts",
    b"survey",
    b"targeted",
    b"tsa",
    b"wgs",
];

/// `s_CompletenessStringToEnum`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/descr_mod_apply.cpp:107-116
const COMPLETENESS: &[&[u8]] = &[
    b"complete",
    b"hasleft",
    b"hasright",
    b"noends",
    b"noleft",
    b"noright",
    b"partial",
];

/// The protein qualifiers of `CFeatModApply::Apply`.
const PROTEIN_QUALS: &[&[u8]] = &[b"protein-desc", b"protein", b"ec-number", b"activity"];

/// C `isspace`.
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// `NStr::IsBlank`.
fn is_blank(text: &[u8]) -> bool {
    text.iter().all(|&c| is_space(c))
}

/// `NStr::Split(.., delims, fSplit_Tokenize)`: the non-empty runs between delimiters.
fn tokenize<'a>(text: &'a [u8], delims: &[u8]) -> Vec<&'a [u8]> {
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
fn normalized_value(value: &[u8]) -> Vec<u8> {
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

/// `CTitleParser::x_FindBrackets`: from `start`, the first `[` that is not inside brackets
/// already open, and its closing `]` (`stop`), with the first `=` inside. `None` where there is
/// no closing bracket, a `]` comes before any `[`, or no `=` comes before the closing bracket.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_reader.cpp:788-824
fn find_brackets(line: &[u8], start: &mut usize) -> Option<(usize, usize)> {
    let mut eq_pos = None;
    let mut open = 0usize;
    let mut i = *start;
    while i < line.len() {
        match line[i] {
            b'[' => {
                open += 1;
                if open == 1 {
                    *start = i;
                }
            }
            b'=' => {
                if open > 0 && eq_pos.is_none() {
                    eq_pos = Some(i);
                }
            }
            b']' => {
                if open == 1 {
                    return eq_pos.filter(|&eq| eq < i).map(|eq| (eq, i));
                } else if open == 0 {
                    return None;
                }
                open -= 1;
            }
            _ => {}
        }
        i += 1;
    }
    None
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
        match find_brackets(title, &mut lb_pos) {
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
    fn single_value<'a>(values: &'a [ModData]) -> &'a [u8] {
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
        match canonical {
            b"location" => {
                if !GENOME.contains(&normalized_value(value).as_slice()) {
                    return fail(&values[0], b"");
                }
                self.effects.non_title = true;
            }
            b"origin" => {
                if !ORIGIN.contains(&normalized_value(value).as_slice()) {
                    return fail(&values[0], b"");
                }
                self.effects.non_title = true;
            }
            b"focus" => {
                if value.eq_ignore_ascii_case(b"true") {
                    self.effects.non_title = true;
                } else if value.eq_ignore_ascii_case(b"false") {
                    return fail(&values[0], b"");
                }
            }
            b"fwd-primer-name" | b"fwd-primer-seq" => self.effects.non_title = true,
            b"rev-primer-name" | b"rev-primer-seq" => {
                if values
                    .iter()
                    .any(|modifier| !tokenize(&modifier.value, b",").is_empty())
                {
                    self.effects.non_title = true;
                }
            }
            b"taxname" | b"common" | b"dbxref" | b"lineage" | b"division" => {
                self.effects.non_title = true
            }
            b"taxid" => {
                if string_to_int(value, 0).is_err() {
                    return fail(&values[0], b"Integer value expected.");
                }
                self.effects.non_title = true;
            }
            b"gcode" | b"mgcode" | b"pgcode" => {
                if string_to_int(value, 0).is_err() {
                    return fail(&values[0], b"Integer value expected.");
                }
                self.effects.non_title = true;
            }
            b"sra" | b"bioproject" | b"biosample" | b"keyword" => {
                if values
                    .iter()
                    .any(|modifier| !tokenize(&modifier.value, b",; \t").is_empty())
                {
                    self.effects.non_title = true;
                }
            }
            b"mol-type" => {
                let Some((_, mol)) = BIOMOL
                    .iter()
                    .find(|(key, _)| *key == normalized_value(value).as_slice())
                else {
                    return fail(&values[0], b"");
                };
                self.effects.non_title = true;
                self.mol = *mol;
                self.mol_set = true;
            }
            b"completeness" => {
                if !COMPLETENESS.contains(&normalized_value(value).as_slice()) {
                    return fail(&values[0], b"");
                }
                self.effects.non_title = true;
            }
            b"tech" => {
                if !TECH.contains(&normalized_value(value).as_slice()) {
                    return fail(&values[0], b"");
                }
                self.effects.non_title = true;
            }
            b"primary-accession" => {
                for modifier in values {
                    if std::str::from_utf8(&modifier.value).is_err() {
                        return Err(ModsStop::Undefined(INVALID_UTF8));
                    }
                }
                if values
                    .iter()
                    .any(|modifier| !tokenize(&modifier.value, b",; \t").is_empty())
                {
                    self.effects.non_title = true;
                }
            }
            b"secondary-accession" => return Err(ModsStop::Undefined(SECONDARY_ACCESSION)),
            b"project" => {
                let mut any = false;
                for modifier in values {
                    let tokens = tokenize(&modifier.value, b",; \t");
                    if tokens.iter().any(|token| string_to_uint(token).is_err()) {
                        return fail(modifier, b"");
                    }
                    any |= !tokens.is_empty();
                }
                if any {
                    self.effects.non_title = true;
                }
            }
            b"comment" | b"ft-map" | b"ft-mod" => self.effects.non_title = true,
            b"pmid" => {
                for modifier in values {
                    if string_to_int8(&modifier.value, 0).is_err() {
                        return fail(&values[0], b"Expected integer value.");
                    }
                }
                self.effects.non_title = true;
            }
            name if SUBSOURCE.contains(&name) => {
                if NEEDS_NO_TEXT.contains(&name) {
                    for modifier in values {
                        if !modifier.value.eq_ignore_ascii_case(b"true") {
                            return fail(modifier, b"");
                        }
                    }
                }
                self.effects.non_title = true;
            }
            name if ORGMOD.contains(&name) => self.effects.non_title = true,
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
                if ![&b"single"[..], b"double", b"mixed", b"other"].contains(&normalized.as_slice())
                {
                    return fail();
                }
                self.effects.non_title = true;
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
                if ![&b"linear"[..], b"circular", b"tandem", b"other"]
                    .contains(&normalized.as_slice())
                {
                    return fail();
                }
                self.effects.non_title = true;
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
            self.effects.non_title = true;
        }
        true
    }
}

/// The cause of the refusal of a secondary accession (`x_SetHist` makes a Seq-id of each value
/// and an invalid one throws a `CSeqIdException` that no `catch` of the reader takes).
pub(crate) const SECONDARY_ACCESSION: &str =
    "NCBI BLAST+ 2.17.0 turns the values of a secondary-accession modifier into Seq-ids, and what it does with one that is not an accession is not followed by LOSAT";

/// The cause of the refusal of a primary accession that is not UTF-8.
pub(crate) const INVALID_UTF8: &str =
    "NCBI BLAST+ 2.17.0 converts the values of a primary-accession modifier to UTF-8, and what it does with one that is not UTF-8 is not followed by LOSAT";

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
