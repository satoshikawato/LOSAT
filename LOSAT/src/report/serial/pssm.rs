//! NCBI-ScoreMat's `PssmWithParameters` (c++/src/objects/scoremat/scoremat.asn): the value that psiblast
//! writes with `-out_pssm` (`MSerial_AsnText`), and its writer over an `Encoder`, in the member
//! order of the ASN.1 module, so that the ASN.1 binary and the archive can write the same value.
//! The reader of the ASN.1 text is `blastinput::pssm_text`.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:385-396
//!
//! A member is written when it is set: an OPTIONAL member that is `None` is not, and neither is a
//! member with a DEFAULT that was not set (`is_protein` and `by_row` hold whether and how they were
//! set). The reader keeps every member NCBI's object keeps (`blastinput::pssm_text`, line RG item
//! D3-1), and the writer writes them all: the members that the search does not read (the
//! descriptors of the query other than its titles, its `annot`, the members of its `Seq-inst` other
//! than `repr`, `mol`, `length` and `seq-data`, a query Seq-entry `set`, Seq-data in a coding
//! other than the six that `SeqData` types, and the `constraints` and `annotatedSites` of the
//! parameters) are values of the
//! generic reader (`report::serial::asn_types::Value`), written by the generic writer
//! (`asn_types::write`, line RG item D3-2).

use super::asn_types::table::named;
use super::asn_types::write::write_value as write_generic;
use super::asn_types::Value;
use super::encoder::Encoder;

/// `PssmWithParameters`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:345-366
#[derive(Debug, Clone, PartialEq)]
pub struct PssmWithParameters {
    pub pssm: Pssm,
    pub params: Option<PssmParameters>,
}

/// `Pssm`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:259-297
#[derive(Debug, Clone, PartialEq)]
pub struct Pssm {
    /// `isProtein`: `Some` when set (written), `None` for the default TRUE (not written).
    pub is_protein: Option<bool>,
    pub identifier: Option<ObjectId>,
    pub num_rows: i32,
    pub num_columns: i32,
    pub row_labels: Option<Vec<Vec<u8>>>,
    /// `byRow`: `Some` when set (written), `None` for the default FALSE (not written).
    pub by_row: Option<bool>,
    pub query: Option<SeqEntry>,
    pub intermediate_data: Option<PssmIntermediateData>,
    pub final_data: Option<PssmFinalData>,
}

/// `Object-id`: the Seq-id model's (`seqid::model`), which the PSSM's identifier shares.
pub use crate::seqid::model::ObjectId;

/// The PSSM's query: the `Seq-entry` `seq` variant, a Bioseq with its ids, its descriptors, its
/// `Seq-inst` and its `annot`.
// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:36-40
#[derive(Debug, Clone, PartialEq)]
pub struct PssmQuery {
    /// The ids, read by the table-driven reader's Seq-id model (`asn_types::seq_id`) and written
    /// as every Seq-id of NCBI's ASN.1 text (`report::seqid::SeqId::Db`).
    pub ids: Vec<crate::seqid::ParsedSeqId>,
    /// `descr` (a `Seq-descr`, a set of `Seqdesc`) in the order read; `None` without the member.
    pub descr: Option<Vec<Seqdesc>>,
    pub inst: SeqInst,
    /// `annot SET OF Seq-annot`, each a value of the generic reader.
    pub annot: Option<Vec<Value>>,
}

/// The PSSM's `query Seq-entry`: its `seq` variant, the Bioseq that the search reads, or its `set`
/// variant, a Bioseq-set as the generic reader read it, which NCBI's object manager takes with its
/// annotations (`pssm_search::annot_index::index_entry`) and NCBI's application then stops at
/// (`CSeq_entry::GetSeq()`, `pssm_search::in_pssm::psi_search_and_object`).
// NCBI reference (598d8ae6): c++/src/objects/seqset/seqset.asn:58-60,63
#[derive(Debug, Clone, PartialEq)]
pub enum SeqEntry {
    Seq(Box<PssmQuery>),
    Set(Value),
}

impl SeqEntry {
    /// The Bioseq of the `seq` variant.
    pub fn seq(&self) -> Option<&PssmQuery> {
        match self {
            SeqEntry::Seq(query) => Some(query.as_ref()),
            SeqEntry::Set(_) => None,
        }
    }

    /// The Bioseq of the `seq` variant, moved out.
    pub fn into_seq(self) -> Option<PssmQuery> {
        match self {
            SeqEntry::Seq(query) => Some(*query),
            SeqEntry::Set(_) => None,
        }
    }

    /// The Bioseq of the `seq` variant, to change.
    pub fn seq_mut(&mut self) -> Option<&mut PssmQuery> {
        match self {
            SeqEntry::Seq(query) => Some(query.as_mut()),
            SeqEntry::Set(_) => None,
        }
    }
}

/// One `Seqdesc` of the query: a `title`, another variant as the generic reader read it (the
/// `Seqdesc` CHOICE value), or an element of the list stated as `NULL` (`Empty`: the empty
/// `CRef<CSeqdesc>` that NCBI's reader leaves for it, `asn_types::Value::Null`).
// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:47-52
#[derive(Debug, Clone, PartialEq)]
pub enum Seqdesc {
    Title(Vec<u8>),
    Other(Value),
    /// An empty reference: NCBI's query set-up dereferences it (`has_empty`), its writers throw
    /// where they would write it (`Encoder::null_element`).
    Empty,
}

impl Seqdesc {
    /// The name of the descriptor's variant (`title`, `molinfo`, ...); `None` for an empty one.
    pub fn variant_name(&self) -> Option<&'static str> {
        match self {
            Seqdesc::Title(_) => Some("title"),
            Seqdesc::Other(Value::Choice(index, _)) => Some(
                super::asn_types::members(super::asn_types::table::named::SEQDESC)
                    [usize::from(*index)]
                .name(),
            ),
            Seqdesc::Other(_) => unreachable!("a Seqdesc is a CHOICE"),
            Seqdesc::Empty => None,
        }
    }
}

impl PssmQuery {
    /// The titles among the descriptors, in order.
    pub fn titles(&self) -> Vec<&[u8]> {
        self.descr
            .iter()
            .flatten()
            .filter_map(|descriptor| match descriptor {
                Seqdesc::Title(title) => Some(title.as_slice()),
                Seqdesc::Other(_) | Seqdesc::Empty => None,
            })
            .collect()
    }

    /// Whether a descriptor is a `molinfo`.
    pub fn has_molinfo(&self) -> bool {
        self.descr
            .iter()
            .flatten()
            .any(|descriptor| descriptor.variant_name() == Some("molinfo"))
    }

    /// Whether a descriptor is empty (`Seqdesc::Empty`).
    pub fn has_empty(&self) -> bool {
        self.descr
            .iter()
            .flatten()
            .any(|descriptor| *descriptor == Seqdesc::Empty)
    }
}

/// `Seq-inst`.
// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:262-297
#[derive(Debug, Clone, PartialEq)]
pub struct SeqInst {
    /// `repr`: the value's name and number (`raw`, 2).
    pub repr: (&'static str, i32),
    /// `mol`: the value's name and number (`aa`, 3).
    pub mol: (&'static str, i32),
    /// `length` (`TSeqPos`, seq.def: read as Uint4).
    pub length: Option<u32>,
    pub fuzz: Option<Value>,
    /// `topology` (DEFAULT linear): the value's number when set.
    pub topology: Option<i32>,
    /// `strand`: the value's number.
    pub strand: Option<i32>,
    pub seq_data: Option<SeqData>,
    pub ext: Option<Value>,
    pub hist: Option<Value>,
}

impl SeqInst {
    /// A `Seq-inst` of `repr`, `mol`, `length` and `seq_data` only.
    pub fn new(
        repr: (&'static str, i32),
        mol: (&'static str, i32),
        length: Option<u32>,
        seq_data: Option<SeqData>,
    ) -> Self {
        SeqInst {
            repr,
            mol,
            length,
            fuzz: None,
            topology: None,
            strand: None,
            seq_data,
            ext: None,
            hist: None,
        }
    }
}

/// The values of `Seq-inst.repr`.
pub const SEQ_INST_REPR: [(&str, i32); 10] = [
    ("not-set", 0),
    ("virtual", 1),
    ("raw", 2),
    ("seg", 3),
    ("const", 4),
    ("ref", 5),
    ("consen", 6),
    ("map", 7),
    ("delta", 8),
    ("other", 255),
];

/// The values of `Seq-inst.mol`.
pub const SEQ_INST_MOL: [(&str, i32); 6] = [
    ("not-set", 0),
    ("dna", 1),
    ("rna", 2),
    ("aa", 3),
    ("na", 4),
    ("other", 255),
];

/// The `Seq-data` variants that the reader types (the three protein codings that the searches
/// read, and the three nucleotide codings that makeprofiledb's database writer converts), and
/// the others as it read them.
// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:346-358
// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:400-405 (`IUPACna` a StringStore,
// `NCBI2na` and `NCBI4na` OCTET STRINGs)
#[derive(Debug, Clone, PartialEq)]
pub enum SeqData {
    /// `iupacna`: IUPACna letters, a StringStore.
    Iupacna(Vec<u8>),
    /// `iupacaa`: IUPACaa letters, a StringStore.
    Iupacaa(Vec<u8>),
    /// `ncbi2na`: NCBI2na bytes, an OCTET STRING.
    Ncbi2na(Vec<u8>),
    /// `ncbi4na`: NCBI4na bytes, an OCTET STRING.
    Ncbi4na(Vec<u8>),
    /// `ncbieaa`: NCBIeaa letters, a StringStore.
    Ncbieaa(Vec<u8>),
    /// `ncbistdaa`: NCBIstdaa codes, an OCTET STRING.
    Ncbistdaa(Vec<u8>),
    /// Another variant, the `Seq-data` CHOICE value of the generic reader: NCBI's object manager
    /// does not read it, LOSAT's search refuses it (`in_pssm::psi_search_and_object`).
    Other(Value),
}

impl SeqData {
    /// The variant's name and its number in the `Seq-data` CHOICE (`iupacna` 0, ...), and the
    /// name of its type; the number of `Other` is the one it was read with.
    pub fn variant(&self) -> (&'static str, u8, &'static str) {
        match self {
            SeqData::Iupacna(_) => ("iupacna", 0, "IUPACna"),
            SeqData::Iupacaa(_) => ("iupacaa", 1, "IUPACaa"),
            SeqData::Ncbi2na(_) => ("ncbi2na", 2, "NCBI2na"),
            SeqData::Ncbi4na(_) => ("ncbi4na", 3, "NCBI4na"),
            SeqData::Ncbieaa(_) => ("ncbieaa", 7, "NCBIeaa"),
            SeqData::Ncbistdaa(_) => ("ncbistdaa", 9, "NCBIstdaa"),
            SeqData::Other(Value::Choice(index, _)) => {
                let member = &super::asn_types::members(named::SEQ_DATA)[usize::from(*index)];
                let tag = u8::try_from(*index).expect("one of the 11 variants");
                (member.name(), tag, "")
            }
            SeqData::Other(_) => unreachable!("a Seq-data is a CHOICE"),
        }
    }
}

/// `PssmIntermediateData`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:197-243
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PssmIntermediateData {
    pub res_freqs_per_pos: Option<Vec<i32>>,
    pub weighted_res_freqs_per_pos: Option<Vec<f64>>,
    pub freq_ratios: Vec<f64>,
    pub information_content: Option<Vec<f64>>,
    pub gapless_column_weights: Option<Vec<f64>>,
    pub sigma: Option<Vec<f64>>,
    pub interval_sizes: Option<Vec<i32>>,
    pub num_matching_seqs: Option<Vec<i32>>,
    pub num_indept_obsr: Option<Vec<f64>>,
}

/// `PssmFinalData`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:154-191
#[derive(Debug, Clone, PartialEq)]
pub struct PssmFinalData {
    pub scores: Vec<i32>,
    pub lambda: f64,
    pub kappa: f64,
    pub h: f64,
    /// `scalingFactor`: `Some` when set, `None` for the default 1 (not written).
    pub scaling_factor: Option<i32>,
    pub lambda_ungapped: Option<f64>,
    pub kappa_ungapped: Option<f64>,
    pub h_ungapped: Option<f64>,
    pub word_score_threshold: Option<f64>,
}

/// `PssmParameters`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:319-341
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PssmParameters {
    pub pseudocount: Option<i32>,
    pub rpsdbparams: Option<FormatRpsDbParameters>,
    /// `constraints CoreDef`, a value of the generic reader.
    pub constraints: Option<Value>,
    pub bit_score_thresh: Option<f64>,
    pub bit_score_reporting_thresh: Option<f64>,
    /// `annotatedSites Site-annot-set`, a value of the generic reader.
    pub annotated_sites: Option<Value>,
}

/// `FormatRpsDbParameters`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:303-315
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FormatRpsDbParameters {
    pub matrix_name: Vec<u8>,
    pub gap_open: Option<i32>,
    pub gap_extend: Option<i32>,
}

/// The six Karlin-Altschul values of a `finalData`, written in place of the object's: psiblast's
/// K/A write-back into the PSSM it searched (`psiblast_impl.cpp:150-168`, cited at
/// `algorithm::psiblast::in_pssm_files::InPssmRounds::written_back`), which sets all six.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FinalKarlin {
    pub lambda: f64,
    pub kappa: f64,
    pub h: f64,
    pub lambda_ungapped: f64,
    pub kappa_ungapped: f64,
    pub h_ungapped: f64,
}

impl PssmWithParameters {
    /// A copy of the object with `karlin` written back into its `finalData` (what the writers'
    /// `karlin` writes in place of the six values; an object without `finalData` stays without):
    /// the PSSM that psiblast's archive of a round holds (`algorithm::blastp::archive::PsiArchive`).
    pub fn with_final_karlin(&self, karlin: FinalKarlin) -> Self {
        let mut value = self.clone();
        if let Some(data) = value.pssm.final_data.as_mut() {
            data.lambda = karlin.lambda;
            data.kappa = karlin.kappa;
            data.h = karlin.h;
            data.lambda_ungapped = Some(karlin.lambda_ungapped);
            data.kappa_ungapped = Some(karlin.kappa_ungapped);
            data.h_ungapped = Some(karlin.h_ungapped);
        }
        value
    }
}

/// The document of a PSSM (`WriteFileHeader`, the object, `EndOfWrite` and the separator), with
/// `karlin` in place of the six values of `finalData`.
fn write_document<E: Encoder>(
    enc: &mut E,
    value: &PssmWithParameters,
    karlin: Option<FinalKarlin>,
) {
    enc.begin_document("PssmWithParameters");
    write_sequence(enc, value, karlin);
    enc.end_document();
}

/// The `PssmWithParameters` SEQUENCE, with `karlin` in place of the six values of `finalData`.
fn write_sequence<E: Encoder>(
    enc: &mut E,
    value: &PssmWithParameters,
    karlin: Option<FinalKarlin>,
) {
    enc.begin_sequence("PssmWithParameters");
    enc.member("pssm", 0);
    write_pssm(enc, &value.pssm, karlin);
    enc.end_member();
    if let Some(params) = &value.params {
        enc.member("params", 1);
        write_params(enc, params);
        enc.end_member();
    }
    enc.end_sequence();
}

// NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:138-142 (the `pssm` variant of the
// request's queries; line RG F1b-2)
/// Writes `value` as the value of a member or variant of an enclosing document, with the
/// enclosing document's encoder: the search strategy's request writes a PSSM query so
/// (`report::archive::writer`, `queries pssm`), and its strings and octets wrap at the depth where
/// they stand, as NCBI's object stream writes them.
pub fn write_value<E: Encoder>(enc: &mut E, value: &PssmWithParameters) {
    write_sequence(enc, value, None);
}

fn write_ints<E: Encoder>(enc: &mut E, values: &[i32]) {
    enc.begin_list(false);
    for &value in values {
        enc.element();
        enc.int(value);
    }
    enc.end_list();
}

fn write_reals<E: Encoder>(enc: &mut E, values: &[f64]) {
    enc.begin_list(false);
    for &value in values {
        enc.element();
        enc.real(value);
    }
    enc.end_list();
}

fn write_pssm<E: Encoder>(enc: &mut E, pssm: &Pssm, karlin: Option<FinalKarlin>) {
    enc.begin_sequence("Pssm");
    if let Some(is_protein) = pssm.is_protein {
        enc.member("isProtein", 0);
        enc.boolean(is_protein);
        enc.end_member();
    }
    if let Some(identifier) = &pssm.identifier {
        enc.member("identifier", 1);
        crate::report::seqid::write_object_id(enc, identifier);
        enc.end_member();
    }
    enc.member("numRows", 2);
    enc.int(pssm.num_rows);
    enc.end_member();
    enc.member("numColumns", 3);
    enc.int(pssm.num_columns);
    enc.end_member();
    if let Some(labels) = &pssm.row_labels {
        enc.member("rowLabels", 4);
        enc.begin_list(false);
        for label in labels {
            enc.element();
            enc.visible_str(label);
        }
        enc.end_list();
        enc.end_member();
    }
    if let Some(by_row) = pssm.by_row {
        enc.member("byRow", 5);
        enc.boolean(by_row);
        enc.end_member();
    }
    if let Some(query) = &pssm.query {
        enc.member("query", 6);
        match query {
            SeqEntry::Seq(query) => {
                enc.variant("Seq-entry", "seq", 0);
                write_query(enc, query);
            }
            // NCBI stops a search of such a query (`in_pssm::psi_search_and_object`) before it
            // writes a PSSM file or a request; the writer writes it as read.
            SeqEntry::Set(set) => {
                enc.variant("Seq-entry", "set", 1);
                write_generic(enc, named::BIOSEQ_SET, set);
            }
        }
        enc.end_variant();
        enc.end_member();
    }
    if let Some(data) = &pssm.intermediate_data {
        enc.member("intermediateData", 7);
        write_intermediate_data(enc, data);
        enc.end_member();
    }
    if let Some(data) = &pssm.final_data {
        enc.member("finalData", 8);
        write_final_data(enc, data, karlin);
        enc.end_member();
    }
    enc.end_sequence();
}

/// The query's Bioseq (`CPsiBlastInputData::x_ExtractQueryForPssm` makes psiblast's).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:144-175
fn write_query<E: Encoder>(enc: &mut E, query: &PssmQuery) {
    enc.begin_sequence("Bioseq");
    enc.member("id", 0);
    enc.begin_list(true);
    for id in &query.ids {
        enc.element();
        crate::report::seqid::SeqId::Db(id).write(enc);
    }
    enc.end_list();
    enc.end_member();
    if let Some(descr) = &query.descr {
        enc.member("descr", 1);
        enc.begin_named("Seq-descr");
        enc.begin_list(true);
        for descriptor in descr {
            match descriptor {
                Seqdesc::Title(title) => {
                    enc.element();
                    enc.variant("Seqdesc", "title", 4);
                    enc.visible_str(title);
                    enc.end_variant();
                }
                Seqdesc::Other(value) => {
                    enc.element();
                    write_generic(enc, named::SEQDESC, value);
                }
                // NCBI's throw (`encoder::NullElement`), before the element's separator.
                Seqdesc::Empty => enc.null_element("Seqdesc"),
            }
        }
        enc.end_list();
        enc.end_named();
        enc.end_member();
    }
    enc.member("inst", 2);
    enc.begin_sequence("Seq-inst");
    enc.member("repr", 0);
    enc.enumerated(query.inst.repr.0, query.inst.repr.1);
    enc.end_member();
    enc.member("mol", 1);
    enc.enumerated(query.inst.mol.0, query.inst.mol.1);
    enc.end_member();
    if let Some(length) = query.inst.length {
        enc.member("length", 2);
        // `WriteUint4`: the decimal digits (and in BER the INTEGER of the value).
        match i32::try_from(length) {
            Ok(length) => enc.int(length),
            Err(_) => enc.big_int(i64::from(length)),
        }
        enc.end_member();
    }
    if let Some(fuzz) = &query.inst.fuzz {
        enc.member("fuzz", 3);
        write_generic(enc, named::INT_FUZZ, fuzz);
        enc.end_member();
    }
    for (name, tag, value) in [
        ("topology", 4, query.inst.topology),
        ("strand", 5, query.inst.strand),
    ] {
        if let Some(value) = value {
            enc.member(name, tag);
            write_generic(
                enc,
                super::asn_types::members(named::SEQ_INST)[usize::from(tag)].ty,
                &Value::Enum(value),
            );
            enc.end_member();
        }
    }
    if let Some(data) = &query.inst.seq_data {
        enc.member("seq-data", 6);
        match data {
            // NCBI stops a search of such a query (`in_pssm::psi_search_and_object`) before it
            // writes a PSSM file or a request; the writer writes it as read. (So does it for a
            // nucleotide coding, which only makeprofiledb reads.)
            SeqData::Other(value) => write_generic(enc, named::SEQ_DATA, value),
            SeqData::Iupacna(letters)
            | SeqData::Iupacaa(letters)
            | SeqData::Ncbieaa(letters)
            | SeqData::Ncbi2na(letters)
            | SeqData::Ncbi4na(letters)
            | SeqData::Ncbistdaa(letters) => {
                let (variant, tag, named_type) = data.variant();
                enc.variant("Seq-data", variant, tag);
                enc.begin_named(named_type);
                match data {
                    SeqData::Iupacna(_) | SeqData::Iupacaa(_) | SeqData::Ncbieaa(_) => {
                        enc.string_store(letters)
                    }
                    _ => enc.octets(letters),
                }
                enc.end_named();
                enc.end_variant();
            }
        }
        enc.end_member();
    }
    for (name, tag, value) in [
        ("ext", 7, (named::SEQ_EXT, &query.inst.ext)),
        ("hist", 8, (named::SEQ_HIST, &query.inst.hist)),
    ] {
        if let (ty, Some(value)) = value {
            enc.member(name, tag);
            write_generic(enc, ty, value);
            enc.end_member();
        }
    }
    enc.end_sequence();
    enc.end_member();
    if let Some(annots) = &query.annot {
        enc.member("annot", 3);
        enc.begin_list(true);
        for annot in annots {
            enc.element();
            write_generic(enc, named::SEQ_ANNOT, annot);
        }
        enc.end_list();
        enc.end_member();
    }
    enc.end_sequence();
}

fn write_intermediate_data<E: Encoder>(enc: &mut E, data: &PssmIntermediateData) {
    enc.begin_sequence("PssmIntermediateData");
    if let Some(values) = &data.res_freqs_per_pos {
        enc.member("resFreqsPerPos", 0);
        write_ints(enc, values);
        enc.end_member();
    }
    if let Some(values) = &data.weighted_res_freqs_per_pos {
        enc.member("weightedResFreqsPerPos", 1);
        write_reals(enc, values);
        enc.end_member();
    }
    enc.member("freqRatios", 2);
    write_reals(enc, &data.freq_ratios);
    enc.end_member();
    for (name, tag, values) in [
        ("informationContent", 3, &data.information_content),
        ("gaplessColumnWeights", 4, &data.gapless_column_weights),
        ("sigma", 5, &data.sigma),
    ] {
        if let Some(values) = values {
            enc.member(name, tag);
            write_reals(enc, values);
            enc.end_member();
        }
    }
    for (name, tag, values) in [
        ("intervalSizes", 6, &data.interval_sizes),
        ("numMatchingSeqs", 7, &data.num_matching_seqs),
    ] {
        if let Some(values) = values {
            enc.member(name, tag);
            write_ints(enc, values);
            enc.end_member();
        }
    }
    if let Some(values) = &data.num_indept_obsr {
        enc.member("numIndeptObsr", 8);
        write_reals(enc, values);
        enc.end_member();
    }
    enc.end_sequence();
}

fn write_final_data<E: Encoder>(enc: &mut E, data: &PssmFinalData, karlin: Option<FinalKarlin>) {
    let (lambda, kappa, h, lambda_ungapped, kappa_ungapped, h_ungapped) = match karlin {
        Some(karlin) => (
            karlin.lambda,
            karlin.kappa,
            karlin.h,
            Some(karlin.lambda_ungapped),
            Some(karlin.kappa_ungapped),
            Some(karlin.h_ungapped),
        ),
        None => (
            data.lambda,
            data.kappa,
            data.h,
            data.lambda_ungapped,
            data.kappa_ungapped,
            data.h_ungapped,
        ),
    };
    enc.begin_sequence("PssmFinalData");
    enc.member("scores", 0);
    write_ints(enc, &data.scores);
    enc.end_member();
    for (name, tag, value) in [("lambda", 1, lambda), ("kappa", 2, kappa), ("h", 3, h)] {
        enc.member(name, tag);
        enc.real(value);
        enc.end_member();
    }
    if let Some(factor) = data.scaling_factor {
        enc.member("scalingFactor", 4);
        enc.int(factor);
        enc.end_member();
    }
    for (name, tag, value) in [
        ("lambdaUngapped", 5, lambda_ungapped),
        ("kappaUngapped", 6, kappa_ungapped),
        ("hUngapped", 7, h_ungapped),
        ("wordScoreThreshold", 8, data.word_score_threshold),
    ] {
        if let Some(value) = value {
            enc.member(name, tag);
            enc.real(value);
            enc.end_member();
        }
    }
    enc.end_sequence();
}

fn write_params<E: Encoder>(enc: &mut E, params: &PssmParameters) {
    enc.begin_sequence("PssmParameters");
    if let Some(pseudocount) = params.pseudocount {
        enc.member("pseudocount", 0);
        enc.int(pseudocount);
        enc.end_member();
    }
    if let Some(rps) = &params.rpsdbparams {
        enc.member("rpsdbparams", 1);
        enc.begin_sequence("FormatRpsDbParameters");
        enc.member("matrixName", 0);
        enc.visible_str(&rps.matrix_name);
        enc.end_member();
        if let Some(gap_open) = rps.gap_open {
            enc.member("gapOpen", 1);
            enc.int(gap_open);
            enc.end_member();
        }
        if let Some(gap_extend) = rps.gap_extend {
            enc.member("gapExtend", 2);
            enc.int(gap_extend);
            enc.end_member();
        }
        enc.end_sequence();
        enc.end_member();
    }
    if let Some(constraints) = &params.constraints {
        enc.member("constraints", 2);
        write_generic(enc, named::COREDEF, constraints);
        enc.end_member();
    }
    for (name, tag, value) in [
        ("bitScoreThresh", 3, params.bit_score_thresh),
        (
            "bitScoreReportingThresh",
            4,
            params.bit_score_reporting_thresh,
        ),
    ] {
        if let Some(value) = value {
            enc.member(name, tag);
            enc.real(value);
            enc.end_member();
        }
    }
    if let Some(sites) = &params.annotated_sites {
        enc.member("annotatedSites", 5);
        write_generic(enc, named::SITE_ANNOT_SET, sites);
        enc.end_member();
    }
    enc.end_sequence();
}

/// `value` in ASN.1 text (`MSerial_AsnText`), as one document; `Err` when the encoder refuses a
/// value (a byte outside the visible ASCII range in a VisibleString: `asn_text_posting` writes it
/// as NCBI does).
#[cfg(test)]
pub fn asn_text(value: &PssmWithParameters) -> Result<Vec<u8>, super::encoder::SerialError> {
    asn_text_with_karlin(value, None)
}

/// `asn_text` of a PSSM that psiblast writes to its file: the document with `#` for each byte of a
/// VisibleString that NCBI's ASN.1 text stream replaces (a query title outside 0x20..0x7E), and
/// NCBI's post of each such byte (`encoder::ReplacedChar`, positions from the document's start),
/// which psiblast makes while it writes the file (`SavePssmToFile`, cited at the module). A PSSM
/// without such a byte is written by `asn_text` alone; one with it again through `Traced`.
///
/// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:413-444 (cited at
/// `seqalign::write_query_posting`: `ReplaceVisibleChar` for each byte that is not
/// `GoodVisibleChar`)
pub fn asn_text_posting(
    value: &PssmWithParameters,
) -> Result<(Vec<u8>, Vec<super::encoder::ReplacedChar>), super::encoder::SerialError> {
    asn_text_posting_with_karlin(value, None)
}

/// `asn_text_posting` of `value` with `karlin`, when given, in place of the six values of its
/// `finalData` (`asn_text_with_karlin`), for a PSSM that LOSAT made (an element stated as `NULL`,
/// which only an `-in_pssm` file holds, is the error `SerialError::NullElement` here:
/// `asn_text_written`).
pub fn asn_text_posting_with_karlin(
    value: &PssmWithParameters,
    karlin: Option<FinalKarlin>,
) -> Result<(Vec<u8>, Vec<super::encoder::ReplacedChar>), super::encoder::SerialError> {
    let written = asn_text_written(value, karlin)?;
    match written.thrown {
        Some(_) => Err(super::encoder::SerialError::NullElement),
        None => Ok((written.bytes, written.replaced)),
    }
}

/// What psiblast's write of a PSSM to its file gives (`asn_text_written`).
#[derive(Debug)]
pub struct WrittenPssm {
    /// The file's bytes: the document, or those that NCBI's stream wrote before it threw.
    pub bytes: Vec<u8>,
    /// NCBI's post of each byte that its stream replaced (`encoder::ReplacedChar`).
    pub replaced: Vec<super::encoder::ReplacedChar>,
    /// NCBI's throw at an element stated as `NULL` (an empty descriptor, `Seqdesc::Empty`, or the
    /// reader's `Value::Null` in another list: `encoder::NullElement`), after the bytes and the
    /// posts.
    pub thrown: Option<Box<super::encoder::NullElement>>,
}

/// `asn_text_posting_with_karlin` of a PSSM that may hold an element stated as `NULL`, as psiblast
/// writes the PSSM of `-in_pssm` (`SavePssmToFile`): a document that `asn_text_with_karlin` does
/// not write as it is (a byte that the stream replaces, an element stated as `NULL`) is written
/// again through
/// `Traced`, which gives NCBI's posts and NCBI's throw.
pub fn asn_text_written(
    value: &PssmWithParameters,
    karlin: Option<FinalKarlin>,
) -> Result<WrittenPssm, super::encoder::SerialError> {
    match asn_text_with_karlin(value, karlin) {
        Err(
            super::encoder::SerialError::NonVisibleChar(_)
            | super::encoder::SerialError::NullElement,
        ) => asn_text_traced(value, karlin),
        written => written.map(|bytes| WrittenPssm {
            bytes,
            replaced: Vec::new(),
            thrown: None,
        }),
    }
}

/// `asn_text_written` through `Traced`.
#[cold]
#[inline(never)]
fn asn_text_traced(
    value: &PssmWithParameters,
    karlin: Option<FinalKarlin>,
) -> Result<WrittenPssm, super::encoder::SerialError> {
    let mut out = Vec::new();
    let mut enc = super::encoder::Traced::new(super::asn_text::AsnTextEncoder::new(&mut out), 0);
    write_document(&mut enc, value, karlin);
    let thrown = enc.take_thrown();
    match enc.error() {
        Some(super::encoder::SerialError::NullElement) if thrown.is_some() => {}
        Some(error) => return Err(error),
        None => {}
    }
    let replaced = enc.take_replaced_chars();
    drop(enc);
    Ok(WrittenPssm {
        bytes: out,
        replaced,
        thrown,
    })
}

/// `asn_text` of `value` with `karlin`, when given, in place of the six values of its `finalData`
/// (the object is not copied).
pub fn asn_text_with_karlin(
    value: &PssmWithParameters,
    karlin: Option<FinalKarlin>,
) -> Result<Vec<u8>, super::encoder::SerialError> {
    let mut out = Vec::new();
    let mut enc = super::asn_text::AsnTextEncoder::new(&mut out);
    write_document(&mut enc, value, karlin);
    if let Some(error) = enc.error() {
        return Err(error);
    }
    Ok(out)
}

#[cfg(test)]
#[path = "pssm_tests.rs"]
mod tests;
