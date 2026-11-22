//! The reader of a PSSM in ASN.1 text or ASN.1 binary (`PssmWithParameters`, the file of
//! psiblast's `-out_pssm` and of `-in_pssm`): NCBI's object streams (`CObjectIStreamAsn`,
//! `CObjectIStreamAsnBinary`) for the types of `report::serial::pssm`, whose writer gives the same
//! bytes back for a PSSM in ASN.1 text that psiblast wrote.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1338-1345
//!
//! The whole document is read by the generic reader of NCBI's type information
//! (`report::serial::asn_types`, line RG item D3-1): every member of every type it reaches, with
//! NCBI's grammar, ranges, replacements and messages (`Post`) and NCBI's errors (`ReadError`). The
//! typed model is then made of the value (`pssm_of`): the members LOSAT's search reads, typed, and
//! the others kept as the generic reader's values: also a query `Seq-entry` that is a `set` and
//! `Seq-data` other than `iupacaa`, `ncbieaa` and `ncbistdaa`, which NCBI's object manager takes
//! (or stops at) before NCBI reads them; the PSSM search refuses them after the object manager's
//! step, for a run that goes on (`pssm_search::in_pssm::psi_search_and_object`). Only the
//! readers' own limits are refused here (`PssmReadError::Refused`).

use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::{
    member_index, members, read_ber_document, read_document, resolved, Kind, Post, ReadError, Value,
};
use crate::report::serial::pssm::{
    FormatRpsDbParameters, ObjectId, Pssm, PssmFinalData, PssmIntermediateData, PssmParameters,
    PssmQuery, PssmWithParameters, SeqData, SeqEntry, SeqInst, Seqdesc,
};

/// Why a PSSM file is not read.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum PssmReadError {
    /// NCBI's reader stops (`ReadError`).
    Ncbi(ReadError),
    /// NCBI reads the PSSM and LOSAT refuses it (what it is, for the refusal's text).
    Refused(String),
}

/// What the reader made of a PSSM file, and the messages NCBI's reader posts while it reads it
/// (before its error, when it stops).
pub(crate) struct PssmRead {
    pub(crate) result: Result<PssmWithParameters, PssmReadError>,
    pub(crate) posts: Vec<Post>,
}

/// Reads the `PssmWithParameters` document at the start of `text` (what follows it is not read,
/// as `stream >> object` reads one object), its values nested to any depth.
pub(crate) fn read_pssm(text: &[u8]) -> PssmRead {
    pssm_read(read_document(text, named::PSSMWITHPARAMETERS))
}

/// [`read_pssm`] of the `PssmWithParameters` object in ASN.1 binary at the start of `bytes`.
pub(crate) fn read_pssm_binary(bytes: &[u8]) -> PssmRead {
    pssm_read(read_ber_document(bytes, named::PSSMWITHPARAMETERS))
}

fn pssm_read((value, posts): (Result<Value, ReadError>, Vec<Post>)) -> PssmRead {
    PssmRead {
        result: value
            .map_err(|error| match error {
                // LOSAT's own limit of a reader: a refusal, not NCBI's error.
                ReadError::Unsupported(what) => PssmReadError::Refused(what),
                error => PssmReadError::Ncbi(error),
            })
            .map(pssm_of),
        posts,
    }
}

/// `read_pssm`'s PSSM, for the tests: an error names NCBI's class or LOSAT's refusal.
#[cfg(test)]
pub(crate) fn read_pssm_with_parameters(text: &[u8]) -> Result<PssmWithParameters, PssmReadError> {
    read_pssm(text).result
}

// The typed boundary: the model of the value read, its parts moved out of it. The value has the
// table's shape (the reader read it under the table), so a mandatory member is there and each
// member has its type's form.

/// The member `name` of the SEQUENCE `value` of the type `ty`, taken out, when set.
fn take(value: &mut Value, ty: u16, name: &str) -> Option<Value> {
    value.take_member(member_index(ty, name))
}

fn int(value: Value) -> i64 {
    match value {
        Value::Int(n) => n,
        other => unreachable!("an INTEGER of the table, not {other:?}"),
    }
}

/// An INTEGER of Int4 storage.
fn int4(value: Value) -> i32 {
    i32::try_from(int(value)).expect("an Int4 of the table")
}

fn boolean(value: Value) -> bool {
    match value {
        Value::Bool(value) => value,
        other => unreachable!("a BOOLEAN of the table, not {other:?}"),
    }
}

fn real(value: Value) -> f64 {
    match value {
        Value::Real(value) => value,
        other => unreachable!("a REAL of the table, not {other:?}"),
    }
}

fn string(mut value: Value) -> Vec<u8> {
    match &mut value {
        Value::Str(text) => std::mem::take(text),
        other => unreachable!("a string of the table, not {other:?}"),
    }
}

fn list(mut value: Value) -> Vec<Value> {
    match &mut value {
        Value::List(items) => std::mem::take(items),
        other => unreachable!("a SEQUENCE OF of the table, not {other:?}"),
    }
}

/// A SEQUENCE OF INTEGER of Int4 storage.
fn int4s(mut value: Value) -> Vec<i32> {
    match &mut value {
        Value::Int4s(values) => std::mem::take(values),
        other => unreachable!("a SEQUENCE OF INTEGER of the table, not {other:?}"),
    }
}

/// A SEQUENCE OF REAL.
fn reals(mut value: Value) -> Vec<f64> {
    match &mut value {
        Value::Reals(values) => std::mem::take(values),
        other => unreachable!("a SEQUENCE OF REAL of the table, not {other:?}"),
    }
}

fn choice(mut value: Value) -> (usize, Value) {
    match &mut value {
        Value::Choice(index, value) => (usize::from(*index), value.take()),
        other => unreachable!("a CHOICE of the table, not {other:?}"),
    }
}

/// An ENUMERATED value: its name in the list of the member's type, and its number.
fn enumerated(value: Value, ty: u16, name: &str) -> (&'static str, i32) {
    let number = enum_number(value);
    let member_ty = members(ty)[member_index(ty, name)].ty;
    let Kind::Enumerated { values, .. } = resolved(member_ty).kind else {
        unreachable!("{name} is an ENUMERATED")
    };
    let name = values
        .values()
        .iter()
        .find(|value| value.value == number)
        .map(|value| value.name())
        .expect("an ENUMERATED number has a name");
    (name, number)
}

fn enum_number(value: Value) -> i32 {
    match value {
        Value::Enum(number) => number,
        other => unreachable!("an ENUMERATED of the table, not {other:?}"),
    }
}

fn object_id(value: Value) -> ObjectId {
    match choice(value) {
        (0, value) => ObjectId::Id(int4(value)),
        (_, value) => ObjectId::Str(string(value)),
    }
}

/// The model of a `PssmWithParameters` value, or LOSAT's refusal.
fn pssm_of(mut value: Value) -> PssmWithParameters {
    const T: u16 = named::PSSMWITHPARAMETERS;
    PssmWithParameters {
        pssm: pssm(take(&mut value, T, "pssm").expect("mandatory")),
        params: take(&mut value, T, "params").map(params),
    }
}

fn pssm(mut value: Value) -> Pssm {
    const T: u16 = named::PSSM;
    let mut get = |name: &str| take(&mut value, T, name);
    Pssm {
        is_protein: get("isProtein").map(boolean),
        identifier: get("identifier").map(object_id),
        num_rows: get("numRows").map(int4).expect("mandatory"),
        num_columns: get("numColumns").map(int4).expect("mandatory"),
        row_labels: get("rowLabels").map(|labels| list(labels).into_iter().map(string).collect()),
        by_row: get("byRow").map(boolean),
        query: get("query").map(seq_entry),
        intermediate_data: get("intermediateData").map(intermediate_data),
        final_data: get("finalData").map(final_data),
    }
}

/// The query's `Seq-entry`: its `seq` variant, a Bioseq, or its `set` variant as read.
fn seq_entry(value: Value) -> SeqEntry {
    let (variant, mut bioseq) = choice(value);
    if members(named::SEQ_ENTRY)[variant].name() != "seq" {
        return SeqEntry::Set(bioseq);
    }
    const T: u16 = named::BIOSEQ;
    let mut get = |name: &str| take(&mut bioseq, T, name);
    SeqEntry::Seq(Box::new(PssmQuery {
        ids: list(get("id").expect("mandatory"))
            .into_iter()
            .map(|id| {
                id.into_seq_id()
                    .unwrap_or_else(|| unreachable!("a Seq-id is read into the Seq-id model"))
            })
            .collect(),
        descr: get("descr").map(|descr| {
            list(descr)
                .into_iter()
                .map(|mut descriptor| {
                    if let Value::Choice(variant, title) = &mut descriptor {
                        if members(named::SEQDESC)[usize::from(*variant)].name() == "title" {
                            return Seqdesc::Title(string(title.take()));
                        }
                    }
                    Seqdesc::Other(descriptor)
                })
                .collect()
        }),
        inst: seq_inst(get("inst").expect("mandatory")),
        annot: get("annot").map(list),
    }))
}

fn seq_inst(mut value: Value) -> SeqInst {
    const T: u16 = named::SEQ_INST;
    let mut get = |name: &str| take(&mut value, T, name);
    let seq_data = match get("seq-data") {
        None => None,
        Some(data) => {
            let (variant, mut value) = choice(data);
            let index = u16::try_from(variant).expect("a CHOICE variant");
            Some(match members(named::SEQ_DATA)[variant].name() {
                "iupacaa" => SeqData::Iupacaa(string(value)),
                "ncbieaa" => SeqData::Ncbieaa(string(value)),
                "ncbistdaa" => SeqData::Ncbistdaa(match &mut value {
                    Value::Octets(codes) => std::mem::take(codes),
                    other => unreachable!("NCBIstdaa is an OCTET STRING, not {other:?}"),
                }),
                _ => SeqData::Other(Value::Choice(index, Box::new(value))),
            })
        }
    };
    SeqInst {
        repr: enumerated(get("repr").expect("mandatory"), T, "repr"),
        mol: enumerated(get("mol").expect("mandatory"), T, "mol"),
        length: get("length").map(|length| u32::try_from(int(length)).expect("a Uint4")),
        fuzz: get("fuzz"),
        topology: get("topology").map(enum_number),
        strand: get("strand").map(enum_number),
        seq_data,
        ext: get("ext"),
        hist: get("hist"),
    }
}

fn intermediate_data(mut value: Value) -> PssmIntermediateData {
    const T: u16 = named::PSSMINTERMEDIATEDATA;
    let mut get = |name: &str| take(&mut value, T, name);
    PssmIntermediateData {
        res_freqs_per_pos: get("resFreqsPerPos").map(int4s),
        weighted_res_freqs_per_pos: get("weightedResFreqsPerPos").map(reals),
        freq_ratios: get("freqRatios").map(reals).expect("mandatory"),
        information_content: get("informationContent").map(reals),
        gapless_column_weights: get("gaplessColumnWeights").map(reals),
        sigma: get("sigma").map(reals),
        interval_sizes: get("intervalSizes").map(int4s),
        num_matching_seqs: get("numMatchingSeqs").map(int4s),
        num_indept_obsr: get("numIndeptObsr").map(reals),
    }
}

fn final_data(mut value: Value) -> PssmFinalData {
    const T: u16 = named::PSSMFINALDATA;
    let mut get = |name: &str| take(&mut value, T, name);
    PssmFinalData {
        scores: get("scores").map(int4s).expect("mandatory"),
        lambda: get("lambda").map(real).expect("mandatory"),
        kappa: get("kappa").map(real).expect("mandatory"),
        h: get("h").map(real).expect("mandatory"),
        scaling_factor: get("scalingFactor").map(int4),
        lambda_ungapped: get("lambdaUngapped").map(real),
        kappa_ungapped: get("kappaUngapped").map(real),
        h_ungapped: get("hUngapped").map(real),
        word_score_threshold: get("wordScoreThreshold").map(real),
    }
}

fn params(mut value: Value) -> PssmParameters {
    const T: u16 = named::PSSMPARAMETERS;
    let mut get = |name: &str| take(&mut value, T, name);
    PssmParameters {
        pseudocount: get("pseudocount").map(int4),
        rpsdbparams: get("rpsdbparams").map(|mut rps| {
            const R: u16 = named::FORMATRPSDBPARAMETERS;
            FormatRpsDbParameters {
                matrix_name: take(&mut rps, R, "matrixName")
                    .map(string)
                    .expect("mandatory"),
                gap_open: take(&mut rps, R, "gapOpen").map(int4),
                gap_extend: take(&mut rps, R, "gapExtend").map(int4),
            }
        }),
        constraints: get("constraints"),
        bit_score_thresh: get("bitScoreThresh").map(real),
        bit_score_reporting_thresh: get("bitScoreReportingThresh").map(real),
        annotated_sites: get("annotatedSites"),
    }
}

#[cfg(test)]
#[path = "pssm_text_tests.rs"]
mod tests;
