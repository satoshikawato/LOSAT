//! What the modifiers of a defline put on the Bioseq of `CAlnReader::GetSeqEntry`: the descriptors
//! of `CDescrModApply` (through its `CDescrCache`), the protein feature of `CFeatModApply` and the
//! `Seq-inst` members of `CModAdder`. The model keeps what NCBI's objects keep (a BioSource, a
//! MolInfo, a GB-block, user objects, comments, publications), in the order the modifiers create
//! them, and writes each as the generic `Value` of its ASN.1 type, so that the PSSM file and the
//! search strategy write it as NCBI's stream does.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/readers/descr_mod_apply.cpp:1002-1259
//! (`CDescrCache`: one descriptor of each kind is made on first use and then reused, a comment and
//! a publication are made anew for each value)

use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::{member_index, members, resolved, Kind, TextRun, Value};

/// A database tag's value (`Object-id`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum Tag {
    Id(i64),
    Str(Vec<u8>),
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct OrgName {
    mods: Option<Vec<(&'static str, Vec<u8>)>>,
    lineage: Option<Vec<u8>>,
    div: Option<Vec<u8>>,
    gcode: Option<i64>,
    mgcode: Option<i64>,
    pgcode: Option<i64>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct OrgRef {
    taxname: Option<Vec<u8>>,
    common: Option<Vec<u8>>,
    db: Option<Vec<(Vec<u8>, Tag)>>,
    orgname: Option<OrgName>,
}

/// One primer of a PCR reaction (`CPCRPrimer`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Primer {
    seq: Option<Vec<u8>>,
    name: Option<Vec<u8>>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Reaction {
    forward: Option<Vec<Primer>>,
    reverse: Option<Vec<Primer>>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Source {
    genome: Option<&'static str>,
    origin: Option<&'static str>,
    org: OrgRef,
    subtype: Option<Vec<(&'static str, Vec<u8>)>>,
    is_focus: bool,
    pcr: Option<Vec<Reaction>>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct MolInfo {
    biomol: Option<&'static str>,
    tech: Option<&'static str>,
    completeness: Option<&'static str>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct GbBlock {
    extra_accessions: Option<Vec<TextRun>>,
    keywords: Option<Vec<Vec<u8>>>,
}

/// A `User-field`'s data, as the modifiers make it.
#[derive(Clone, Debug, PartialEq, Eq)]
enum FieldData {
    Strs(Vec<Vec<u8>>),
    Str(Vec<u8>),
    Int(i64),
    Fields(Vec<Field>),
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Field {
    /// `label`: a string, or `id 0`.
    label: Option<Vec<u8>>,
    num: Option<i64>,
    data: FieldData,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum UserKind {
    DbLink,
    Tpa,
    GenomeProjects,
    FileTrack,
}

impl UserKind {
    fn type_name(self) -> &'static [u8] {
        match self {
            UserKind::DbLink => b"DBLink",
            UserKind::Tpa => b"TpaAssembly",
            UserKind::GenomeProjects => b"GenomeProjectsDB",
            UserKind::FileTrack => b"FileTrack",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum Item {
    Source(Box<Source>),
    MolInfo(MolInfo),
    Genbank(GbBlock),
    User(UserKind, Vec<Field>),
    Comment(Vec<u8>),
    Pub(i64),
}

/// The protein feature of `CFeatModApply`: a `Prot-ref` on the whole sequence.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct ProtRef {
    pub(crate) name: Option<Vec<Vec<u8>>>,
    pub(crate) desc: Option<Vec<u8>>,
    pub(crate) ec: Option<Vec<Vec<u8>>>,
    pub(crate) activity: Option<Vec<Vec<u8>>>,
}

/// Everything the modifiers of one defline set on the Bioseq besides its title and its molecule.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Content {
    items: Vec<Item>,
    pub(crate) protein: Option<ProtRef>,
    /// `Seq-inst.strand` (`ss`, `ds`, `mixed`, `other`).
    pub(crate) strand: Option<&'static str>,
    /// `Seq-inst.topology`.
    pub(crate) topology: Option<&'static str>,
    /// The ids that `secondary-accession` sets as `Seq-inst.hist.replaces`.
    pub(crate) replaces: Option<Vec<TextRun>>,
    preserve_taxid: bool,
}

fn visible_str(text: &[u8]) -> Value {
    Value::Str(text.to_vec())
}

fn str_list(items: &[Vec<u8>]) -> Value {
    Value::List(items.iter().map(|item| visible_str(item)).collect())
}

/// The type of the member `name` of the SEQUENCE or CHOICE `ty`.
fn child(ty: u16, name: &str) -> u16 {
    members(ty)[member_index(ty, name)].ty
}

/// A SEQUENCE value of `ty` from its set members, in module order.
fn seq(ty: u16, fields: Vec<(&str, Value)>) -> Value {
    let mut set: Vec<(u16, Value)> = fields
        .into_iter()
        .map(|(name, value)| (member_index(ty, name) as u16, value))
        .collect();
    set.sort_by_key(|(position, _)| *position);
    Value::Sequence(set)
}

/// A CHOICE value of `ty`.
fn choice(ty: u16, name: &str, value: Value) -> Value {
    Value::Choice(member_index(ty, name) as u16, Box::new(value))
}

/// The ENUMERATED value `name` of the member `member` of `ty`. Every name that reaches it was
/// checked against the enumeration when the modifier was applied (`mods::checked`, which ends the
/// run with an error where a table of LOSAT's names a value that NCBI's type does not have), so
/// no run can end here.
fn enum_value(ty: u16, member: &str, name: &str) -> Value {
    let found = match &resolved(child(ty, member)).kind {
        Kind::Enumerated { values, .. } => {
            values.values().iter().find(|value| value.name() == name)
        }
        _ => None,
    };
    match found {
        Some(found) => Value::Enum(found.value),
        None => unreachable!("{member} has no value {name}: checked when the modifier was applied"),
    }
}

/// `Object-id`: `id` or `str`.
fn object_id(ty: u16, tag: &Tag) -> Value {
    match tag {
        Tag::Id(id) => choice(ty, "id", Value::Int(*id)),
        Tag::Str(text) => choice(ty, "str", visible_str(text)),
    }
}

impl OrgName {
    fn value(&self, ty: u16) -> Value {
        let mut fields = Vec::new();
        if let Some(mods) = &self.mods {
            let mod_ty = match &resolved(child(ty, "mod")).kind {
                Kind::List { element, .. } => *element,
                _ => unreachable!("OrgName.mod is a list"),
            };
            fields.push((
                "mod",
                Value::List(
                    mods.iter()
                        .map(|(subtype, subname)| {
                            seq(
                                mod_ty,
                                vec![
                                    ("subtype", enum_value(mod_ty, "subtype", subtype)),
                                    ("subname", visible_str(subname)),
                                ],
                            )
                        })
                        .collect(),
                ),
            ));
        }
        if let Some(lineage) = &self.lineage {
            fields.push(("lineage", visible_str(lineage)));
        }
        if let Some(gcode) = self.gcode {
            fields.push(("gcode", Value::Int(gcode)));
        }
        if let Some(mgcode) = self.mgcode {
            fields.push(("mgcode", Value::Int(mgcode)));
        }
        if let Some(div) = &self.div {
            fields.push(("div", visible_str(div)));
        }
        if let Some(pgcode) = self.pgcode {
            fields.push(("pgcode", Value::Int(pgcode)));
        }
        seq(ty, fields)
    }
}

fn list_element(ty: u16, member: &str) -> u16 {
    match &resolved(child(ty, member)).kind {
        Kind::List { element, .. } => *element,
        _ => panic!("{member} is not a list"),
    }
}

impl OrgRef {
    fn value(&self, ty: u16) -> Value {
        let mut fields = Vec::new();
        if let Some(taxname) = &self.taxname {
            fields.push(("taxname", visible_str(taxname)));
        }
        if let Some(common) = &self.common {
            fields.push(("common", visible_str(common)));
        }
        if let Some(db) = &self.db {
            let dbtag = list_element(ty, "db");
            let tag_ty = child(dbtag, "tag");
            fields.push((
                "db",
                Value::List(
                    db.iter()
                        .map(|(name, tag)| {
                            seq(
                                dbtag,
                                vec![("db", visible_str(name)), ("tag", object_id(tag_ty, tag))],
                            )
                        })
                        .collect(),
                ),
            ));
        }
        if let Some(orgname) = &self.orgname {
            fields.push(("orgname", orgname.value(child(ty, "orgname"))));
        }
        seq(ty, fields)
    }
}

fn primer_set(ty: u16, member: &str, primers: &[Primer]) -> Value {
    let primer_ty = match &resolved(child(ty, member)).kind {
        Kind::List { element, .. } => *element,
        _ => panic!("{member} is not a list"),
    };
    Value::List(
        primers
            .iter()
            .map(|primer| {
                let mut fields = Vec::new();
                if let Some(seq_text) = &primer.seq {
                    fields.push(("seq", visible_str(seq_text)));
                }
                if let Some(name) = &primer.name {
                    fields.push(("name", visible_str(name)));
                }
                seq(primer_ty, fields)
            })
            .collect(),
    )
}

impl Source {
    fn value(&self, ty: u16) -> Value {
        let mut fields = Vec::new();
        if let Some(genome) = self.genome {
            fields.push(("genome", enum_value(ty, "genome", genome)));
        }
        if let Some(origin) = self.origin {
            fields.push(("origin", enum_value(ty, "origin", origin)));
        }
        fields.push(("org", self.org.value(child(ty, "org"))));
        if let Some(subtype) = &self.subtype {
            let sub_ty = list_element(ty, "subtype");
            fields.push((
                "subtype",
                Value::List(
                    subtype
                        .iter()
                        .map(|(kind, name)| {
                            seq(
                                sub_ty,
                                vec![
                                    ("subtype", enum_value(sub_ty, "subtype", kind)),
                                    ("name", visible_str(name)),
                                ],
                            )
                        })
                        .collect(),
                ),
            ));
        }
        if self.is_focus {
            fields.push(("is-focus", Value::Null));
        }
        if let Some(pcr) = &self.pcr {
            let reaction_ty = list_element(ty, "pcr-primers");
            fields.push((
                "pcr-primers",
                Value::List(
                    pcr.iter()
                        .map(|reaction| {
                            let mut parts = Vec::new();
                            if let Some(forward) = &reaction.forward {
                                parts
                                    .push(("forward", primer_set(reaction_ty, "forward", forward)));
                            }
                            if let Some(reverse) = &reaction.reverse {
                                parts
                                    .push(("reverse", primer_set(reaction_ty, "reverse", reverse)));
                            }
                            seq(reaction_ty, parts)
                        })
                        .collect(),
                ),
            ));
        }
        seq(ty, fields)
    }
}

fn user_field(ty: u16, field: &Field) -> Value {
    let label_ty = child(ty, "label");
    let data_ty = child(ty, "data");
    let mut fields = Vec::new();
    fields.push((
        "label",
        match &field.label {
            Some(text) => choice(label_ty, "str", visible_str(text)),
            None => choice(label_ty, "id", Value::Int(0)),
        },
    ));
    if let Some(num) = field.num {
        fields.push(("num", Value::Int(num)));
    }
    fields.push((
        "data",
        match &field.data {
            FieldData::Strs(items) => choice(data_ty, "strs", str_list(items)),
            FieldData::Str(text) => choice(data_ty, "str", visible_str(text)),
            FieldData::Int(number) => choice(data_ty, "int", Value::Int(*number)),
            FieldData::Fields(inner) => choice(
                data_ty,
                "fields",
                Value::List(inner.iter().map(|inner| user_field(ty, inner)).collect()),
            ),
        },
    ));
    seq(ty, fields)
}

impl Item {
    fn value(&self) -> Value {
        let ty = named::SEQDESC;
        match self {
            Item::Source(source) => choice(ty, "source", source.value(child(ty, "source"))),
            Item::MolInfo(info) => {
                let info_ty = child(ty, "molinfo");
                let mut fields = Vec::new();
                if let Some(biomol) = info.biomol {
                    fields.push(("biomol", enum_value(info_ty, "biomol", biomol)));
                }
                if let Some(tech) = info.tech {
                    fields.push(("tech", enum_value(info_ty, "tech", tech)));
                }
                if let Some(completeness) = info.completeness {
                    fields.push((
                        "completeness",
                        enum_value(info_ty, "completeness", completeness),
                    ));
                }
                choice(ty, "molinfo", seq(info_ty, fields))
            }
            Item::Genbank(block) => {
                let block_ty = child(ty, "genbank");
                let mut fields = Vec::new();
                if let Some(extra) = &block.extra_accessions {
                    fields.push(("extra-accessions", Value::Texts(extra.clone())));
                }
                if let Some(keywords) = &block.keywords {
                    fields.push(("keywords", str_list(keywords)));
                }
                choice(ty, "genbank", seq(block_ty, fields))
            }
            Item::User(kind, data) => {
                let user_ty = child(ty, "user");
                let type_ty = child(user_ty, "type");
                let field_ty = list_element(user_ty, "data");
                choice(
                    ty,
                    "user",
                    seq(
                        user_ty,
                        vec![
                            (
                                "type",
                                choice(type_ty, "str", visible_str(kind.type_name())),
                            ),
                            (
                                "data",
                                Value::List(
                                    data.iter()
                                        .map(|field| user_field(field_ty, field))
                                        .collect(),
                                ),
                            ),
                        ],
                    ),
                )
            }
            Item::Comment(text) => choice(ty, "comment", visible_str(text)),
            Item::Pub(pmid) => {
                let pubdesc = child(ty, "pub");
                let pub_ty = list_element(pubdesc, "pub");
                choice(
                    ty,
                    "pub",
                    seq(
                        pubdesc,
                        vec![(
                            "pub",
                            Value::List(vec![choice(pub_ty, "pmid", Value::Int(*pmid))]),
                        )],
                    ),
                )
            }
        }
    }
}

impl Content {
    /// The descriptors, in the order the modifiers made them, as `Seqdesc` values.
    pub(crate) fn descriptors(&self) -> Vec<Value> {
        self.items.iter().map(Item::value).collect()
    }

    #[cfg(test)]
    pub(crate) fn has_any(&self) -> bool {
        !self.items.is_empty()
            || self.protein.is_some()
            || self.strand.is_some()
            || self.topology.is_some()
            || self.replaces.is_some()
    }

    fn slot(&mut self, make: impl Fn(&Item) -> bool, new: Item) -> usize {
        match self.items.iter().position(make) {
            Some(at) => at,
            None => {
                self.items.push(new);
                self.items.len() - 1
            }
        }
    }

    fn source(&mut self) -> &mut Source {
        let at = self.slot(
            |item| matches!(item, Item::Source(_)),
            Item::Source(Box::default()),
        );
        match &mut self.items[at] {
            Item::Source(source) => source,
            _ => unreachable!(),
        }
    }

    fn mol_info(&mut self) -> &mut MolInfo {
        let at = self.slot(
            |item| matches!(item, Item::MolInfo(_)),
            Item::MolInfo(MolInfo::default()),
        );
        match &mut self.items[at] {
            Item::MolInfo(info) => info,
            _ => unreachable!(),
        }
    }

    fn genbank(&mut self) -> &mut GbBlock {
        let at = self.slot(
            |item| matches!(item, Item::Genbank(_)),
            Item::Genbank(GbBlock::default()),
        );
        match &mut self.items[at] {
            Item::Genbank(block) => block,
            _ => unreachable!(),
        }
    }

    fn user(&mut self, kind: UserKind) -> &mut Vec<Field> {
        let at = self.slot(
            |item| matches!(item, Item::User(found, _) if *found == kind),
            Item::User(kind, Vec::new()),
        );
        match &mut self.items[at] {
            Item::User(_, fields) => fields,
            _ => unreachable!(),
        }
    }

    // ---- BioSource ----

    pub(crate) fn set_genome(&mut self, name: &'static str) {
        self.source().genome = Some(name);
    }

    pub(crate) fn set_origin(&mut self, name: &'static str) {
        self.source().origin = Some(name);
    }

    pub(crate) fn set_focus(&mut self) {
        self.source().is_focus = true;
    }

    /// `x_SetSubtype`: the subsource of each value (the list is made on the first call).
    pub(crate) fn add_subsources(&mut self, kind: &'static str, values: &[Vec<u8>]) {
        if kind == "plasmid-name" {
            self.source().genome = Some("plasmid");
        }
        let list = self.source().subtype.get_or_insert_with(Vec::new);
        for value in values {
            list.push((kind, value.clone()));
        }
    }

    /// `x_SetOrgMod`: the orgmod of each value (the list is made on the first call).
    pub(crate) fn add_orgmods(&mut self, kind: &'static str, values: &[Vec<u8>]) {
        let orgname = self
            .source()
            .org
            .orgname
            .get_or_insert_with(OrgName::default);
        let list = orgname.mods.get_or_insert_with(Vec::new);
        for value in values {
            list.push((kind, value.clone()));
        }
    }

    pub(crate) fn set_common(&mut self, value: &[u8]) {
        self.source().org.common = Some(value.to_vec());
    }

    pub(crate) fn set_lineage(&mut self, value: &[u8]) {
        self.source()
            .org
            .orgname
            .get_or_insert_with(OrgName::default)
            .lineage = Some(value.to_vec());
    }

    pub(crate) fn set_division(&mut self, value: &[u8]) {
        self.source()
            .org
            .orgname
            .get_or_insert_with(OrgName::default)
            .div = Some(value.to_vec());
    }

    /// `gcode`, `mgcode` or `pgcode`.
    pub(crate) fn set_code(&mut self, which: &[u8], code: i64) {
        let orgname = self
            .source()
            .org
            .orgname
            .get_or_insert_with(OrgName::default);
        match which {
            b"gcode" => orgname.gcode = Some(code),
            b"mgcode" => orgname.mgcode = Some(code),
            _ => orgname.pgcode = Some(code),
        }
    }

    /// `COrg_ref::GetTaxId` / `SetTaxId`: the `taxon` database tag of the organism.
    fn taxon_id(&mut self) -> Option<i64> {
        self.source()
            .org
            .db
            .as_ref()?
            .iter()
            .find_map(|(name, tag)| match (name.as_slice() == b"taxon", tag) {
                (true, Tag::Id(id)) => Some(*id),
                _ => None,
            })
    }

    fn set_taxon_id(&mut self, id: i64) {
        let db = self.source().org.db.get_or_insert_with(Vec::new);
        if let Some((_, tag)) = db.iter_mut().find(|(name, _)| name.as_slice() == b"taxon") {
            *tag = Tag::Id(id);
        } else {
            db.push((b"taxon".to_vec(), Tag::Id(id)));
        }
    }

    /// `taxname`: the name, and a taxon id that no `taxid` modifier set is cleared to 0.
    pub(crate) fn set_taxname(&mut self, value: &[u8]) {
        self.source().org.taxname = Some(value.to_vec());
        if !self.preserve_taxid && self.taxon_id().is_some_and(|id| id != 0) {
            self.set_taxon_id(0);
        }
    }

    pub(crate) fn set_taxid(&mut self, id: i64) {
        self.set_taxon_id(id);
        self.preserve_taxid = true;
    }

    /// `x_SetDBxref`: the database tags, in place of the organism's.
    pub(crate) fn set_dbxref(&mut self, values: &[Vec<u8>]) {
        let tags = values
            .iter()
            .map(|value| match value.iter().position(|&c| c == b':') {
                Some(colon) if colon + 1 < value.len() => (
                    value[..colon].to_vec(),
                    Tag::Str(value[colon + 1..].to_vec()),
                ),
                _ => (b"?".to_vec(), Tag::Str(value.clone())),
            })
            .collect();
        self.source().org.db = Some(tags);
    }

    // ---- PCR primers ----

    fn reactions(&mut self) -> &mut Vec<Reaction> {
        self.source().pcr.get_or_insert_with(Vec::new)
    }

    /// `fwd-primer-name`, `fwd-primer-seq`, `rev-primer-name` and `rev-primer-seq` over the
    /// comma-separated groups of the values.
    pub(crate) fn set_primers(&mut self, forward: bool, names: bool, values: &[Vec<u8>]) {
        let mut groups: Vec<Vec<u8>> = Vec::new();
        for value in values {
            let mut parts: Vec<Vec<u8>> = super::mods::tokenize(value, b",")
                .into_iter()
                .map(<[u8]>::to_vec)
                .collect();
            if !names && parts.len() > 1 {
                if parts[0].first() == Some(&b'(') {
                    parts[0].remove(0);
                }
                let last = parts.len() - 1;
                if parts[last].last() == Some(&b')') {
                    parts[last].pop();
                }
            }
            for part in parts {
                groups.push(if names {
                    part
                } else {
                    part.to_ascii_lowercase()
                });
            }
        }
        if forward {
            let reactions = self.reactions();
            for (at, group) in groups.iter().enumerate() {
                if at >= reactions.len() {
                    reactions.push(Reaction::default());
                }
                set_primer_texts(&mut reactions[at], true, names, group);
            }
        } else if !groups.is_empty() {
            let reactions = self.reactions();
            let have = reactions.len();
            if groups.len() <= have {
                // the last reaction takes the last group, and so on backwards
                let first = have - groups.len();
                for (at, group) in groups.iter().enumerate() {
                    set_primer_texts(&mut reactions[first + at], false, names, group);
                }
            } else {
                for (at, group) in groups.iter().enumerate() {
                    if at >= reactions.len() {
                        reactions.push(Reaction::default());
                    }
                    set_primer_texts(&mut reactions[at], false, names, group);
                }
            }
        }
    }

    // ---- MolInfo ----

    pub(crate) fn set_biomol(&mut self, name: &'static str) {
        self.mol_info().biomol = Some(name);
    }

    pub(crate) fn set_tech(&mut self, name: &'static str) {
        self.mol_info().tech = Some(name);
    }

    pub(crate) fn set_completeness(&mut self, name: &'static str) {
        self.mol_info().completeness = Some(name);
    }

    // ---- GB-block ----

    pub(crate) fn set_keywords(&mut self, keywords: Vec<Vec<u8>>) {
        self.genbank().keywords = Some(keywords);
    }

    pub(crate) fn set_extra_accessions(&mut self, ids: Vec<TextRun>) {
        self.genbank().extra_accessions = Some(ids);
    }

    // ---- user objects ----

    /// `x_SetDBLinkFieldVals`: the field of `label` (found without regard to case) takes `values`.
    pub(crate) fn set_dblink(&mut self, label: &[u8], values: Vec<Vec<u8>>) {
        let fields = self.user(UserKind::DbLink);
        match fields.iter_mut().find(|field| {
            field
                .label
                .as_ref()
                .is_some_and(|found| found.eq_ignore_ascii_case(label))
        }) {
            Some(field) => field.data = FieldData::Strs(values),
            None => fields.push(Field {
                label: Some(label.to_vec()),
                num: None,
                data: FieldData::Strs(values),
            }),
        }
    }

    pub(crate) fn add_file_track(&mut self, label: &[u8], value: &[u8]) {
        self.user(UserKind::FileTrack).push(Field {
            label: Some(label.to_vec()),
            num: Some(1),
            data: FieldData::Str(value.to_vec()),
        });
    }

    pub(crate) fn set_tpa(&mut self, accessions: Vec<Vec<u8>>) {
        *self.user(UserKind::Tpa) = accessions
            .into_iter()
            .map(|accession| Field {
                label: None,
                num: None,
                data: FieldData::Fields(vec![Field {
                    label: Some(b"accession".to_vec()),
                    num: None,
                    data: FieldData::Str(accession),
                }]),
            })
            .collect();
    }

    pub(crate) fn set_genome_projects(&mut self, ids: Vec<i64>) {
        *self.user(UserKind::GenomeProjects) = ids
            .into_iter()
            .map(|id| Field {
                label: None,
                num: None,
                data: FieldData::Fields(vec![
                    Field {
                        label: Some(b"ProjectID".to_vec()),
                        num: None,
                        data: FieldData::Int(id),
                    },
                    Field {
                        label: Some(b"ParentID".to_vec()),
                        num: None,
                        data: FieldData::Int(0),
                    },
                ]),
            })
            .collect();
    }

    // ---- comments and publications ----

    pub(crate) fn add_comment(&mut self, text: &[u8]) {
        self.items.push(Item::Comment(text.to_vec()));
    }

    pub(crate) fn add_pmid(&mut self, pmid: i64) {
        self.items.push(Item::Pub(pmid));
    }

    // ---- the protein feature ----

    pub(crate) fn protein(&mut self) -> &mut ProtRef {
        self.protein.get_or_insert_with(ProtRef::default)
    }
}

/// `s_SetPrimerNames` / `s_SetPrimerSeqs`: the colon-separated primers of `group` into the
/// reaction's forward or reverse set, a blank one skipped.
fn set_primer_texts(reaction: &mut Reaction, forward: bool, names: bool, group: &[u8]) {
    let set = if forward {
        &mut reaction.forward
    } else {
        &mut reaction.reverse
    }
    .get_or_insert_with(Vec::new);
    let size = set.len();
    let mut at = 0;
    for (i, text) in group
        .split(|&c| c == b':')
        .filter(|t| !t.is_empty())
        .enumerate()
    {
        // a blank text is skipped
        if super::mods::is_blank(text) {
            continue;
        }
        if i < size {
            if names {
                set[at].name = Some(text.to_vec());
            } else {
                set[at].seq = Some(text.to_vec());
            }
            at += 1;
        } else {
            let mut primer = Primer::default();
            if names {
                primer.name = Some(text.to_vec());
            } else {
                primer.seq = Some(text.to_vec());
            }
            set.push(primer);
        }
    }
}

/// The ENUMERATED value of the member `member` of `BioSource` whose `g_GetNormalizedModVal` form
/// is `key`, or one of `extras` (NCBI's `s_InitModStringToEnumMap`: the type's own values first).
pub(crate) fn biosource_enum(
    member: &str,
    key: &[u8],
    extras: &[(&[u8], &'static str)],
) -> Option<&'static str> {
    let Kind::Enumerated { values, .. } = &resolved(child(named::BIOSOURCE, member)).kind else {
        return None;
    };
    values
        .values()
        .iter()
        .map(|value| value.name())
        .find(|name| super::mods::normalized_value(name.as_bytes()) == key)
        .or_else(|| {
            extras
                .iter()
                .find(|(extra, _)| super::mods::normalized_value(extra) == key)
                .map(|(_, name)| *name)
        })
}

/// The name of the value `name` of the enumeration `member` of the type `ty`, as the table spells
/// it (`None` for a name that is not one of its values).
pub(crate) fn enum_name(ty: u16, member: &str, name: &[u8]) -> Option<&'static str> {
    let Kind::Enumerated { values, .. } = &resolved(child(ty, member)).kind else {
        return None;
    };
    values
        .values()
        .iter()
        .map(|value| value.name())
        .find(|found| found.as_bytes() == name)
}

/// The type of `SubSource` and `OrgMod`, for [`enum_name`].
pub(crate) fn subsource_type() -> u16 {
    list_element(named::BIOSOURCE, "subtype")
}

pub(crate) fn orgmod_type() -> u16 {
    list_element(child(child(named::BIOSOURCE, "org"), "orgname"), "mod")
}

impl Content {
    /// `Bioseq.annot` of the protein feature (`x_CreateSeqfeat`): one `Seq-annot` with the feature
    /// on the whole sequence, located by the best id of the Bioseq (`FindBestChoice`,
    /// `CSeq_id::BestRank`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/feature_mod_apply.cpp:158-211
    pub(crate) fn annot(&self, ids: &[crate::seqid::ParsedSeqId]) -> Option<Vec<Value>> {
        let prot = self.protein.as_ref()?;
        let annot_ty = named::SEQ_ANNOT;
        let feat_ty = list_element(child(annot_ty, "data"), "ftable");
        let data_ty = child(feat_ty, "data");
        let prot_ty = child(data_ty, "prot");
        let loc_ty = child(feat_ty, "location");
        let mut parts = Vec::new();
        if let Some(name) = &prot.name {
            parts.push(("name", str_list(name)));
        }
        if let Some(desc) = &prot.desc {
            parts.push(("desc", visible_str(desc)));
        }
        if let Some(ec) = &prot.ec {
            parts.push(("ec", str_list(ec)));
        }
        if let Some(activity) = &prot.activity {
            parts.push(("activity", str_list(activity)));
        }
        let best = crate::seqid::rank::find_best(ids, crate::seqid::rank::best_rank_score);
        let mut feat = vec![("data", choice(data_ty, "prot", seq(prot_ty, parts)))];
        if let Some(best) = best {
            feat.push((
                "location",
                choice(loc_ty, "whole", Value::SeqId(Box::new(best.clone()))),
            ));
        }
        let feature = seq(feat_ty, feat);
        Some(vec![seq(
            annot_ty,
            vec![(
                "data",
                choice(
                    child(annot_ty, "data"),
                    "ftable",
                    Value::List(vec![feature]),
                ),
            )],
        )])
    }
}

impl Content {
    /// `Seq-inst.hist` of the `secondary-accession` ids (`x_SetHist`): `replaces` with the ids.
    pub(crate) fn hist(&self) -> Option<Value> {
        let ids = self.replaces.as_ref()?;
        let hist_ty = child(named::SEQ_INST, "hist");
        let rec_ty = child(hist_ty, "replaces");
        Some(seq(
            hist_ty,
            vec![(
                "replaces",
                seq(rec_ty, vec![("ids", Value::Texts(ids.clone()))]),
            )],
        ))
    }
}
