//! The reader of a PSSM in ASN.1 text (`PssmWithParameters`, the file of psiblast's `-out_pssm`
//! and of `-in_pssm`): NCBI's ASN.1 text object stream (`CObjectIStreamAsn`) for the types of
//! `report::serial::pssm`, whose writer gives the same bytes back for a PSSM that psiblast wrote.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1209-1216
//! ```c++
//!     if (args.Exist(kArgPSIInputChkPntFile) && args[kArgPSIInputChkPntFile]) {
//!         CNcbiIstream& in = args[kArgPSIInputChkPntFile].AsInputFile();
//!         _ASSERT(m_Pssm.Empty());
//!         m_Pssm.Reset(new CPssmWithParameters);
//!         try {
//!             switch (CFormatGuess().Format(in)) {
//!             case CFormatGuess::eBinaryASN:
//!                 in >> MSerial_AsnBinary >> *m_Pssm;
//! ```
//!
//! The whole document is read by the generic reader of NCBI's type information
//! (`report::serial::asn_types`, line RG item D3-1): every member of every type it reaches, with
//! NCBI's grammar, ranges, replacements and messages (`Post`) and NCBI's errors (`ReadError`). The
//! typed model is then made of the value (`pssm_of`): the members LOSAT's search reads, typed, and
//! the others kept as the generic reader's values. LOSAT refuses at that boundary, after NCBI's
//! reader has read everything, what its PSSM search does not take (`PssmReadError::Refused`): a
//! query `Seq-entry` that is a `set`, `Seq-data` other than `iupacaa`, `ncbieaa` and `ncbistdaa`.

use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::{
    member_index, members, read_document, resolved, Kind, Post, ReadError, Value,
};
use crate::report::serial::pssm::{
    FormatRpsDbParameters, ObjectId, Pssm, PssmFinalData, PssmIntermediateData, PssmParameters,
    PssmQuery, PssmWithParameters, SeqData, SeqInst, Seqdesc,
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
/// as `stream >> object` reads one object).
pub(crate) fn read_pssm(text: &[u8]) -> PssmRead {
    let (value, posts) = read_document(text, named::PSSMWITHPARAMETERS);
    PssmRead {
        result: value
            .map_err(PssmReadError::Ncbi)
            .and_then(|value| pssm_of(&value).map_err(PssmReadError::Refused)),
        posts,
    }
}

/// `read_pssm`'s PSSM, for the tests: an error names NCBI's class or LOSAT's refusal.
#[cfg(test)]
pub(crate) fn read_pssm_with_parameters(text: &[u8]) -> Result<PssmWithParameters, PssmReadError> {
    read_pssm(text).result
}

// The typed boundary: the model of the value read. The value has the table's shape (the reader
// read it under the table), so a mandatory member is there and each member has its type's form.

fn member<'v>(value: &'v Value, ty: u16, name: &str) -> Option<&'v Value> {
    value.member(member_index(ty, name))
}

fn int(value: &Value) -> i64 {
    match value {
        Value::Int(n) => *n,
        other => unreachable!("an INTEGER of the table, not {other:?}"),
    }
}

/// An INTEGER of Int4 storage.
fn int4(value: &Value) -> i32 {
    i32::try_from(int(value)).expect("an Int4 of the table")
}

fn boolean(value: &Value) -> bool {
    match value {
        Value::Bool(value) => *value,
        other => unreachable!("a BOOLEAN of the table, not {other:?}"),
    }
}

fn real(value: &Value) -> f64 {
    match value {
        Value::Real(value) => *value,
        other => unreachable!("a REAL of the table, not {other:?}"),
    }
}

fn string(value: &Value) -> Vec<u8> {
    match value {
        Value::Str(text) => text.clone(),
        other => unreachable!("a string of the table, not {other:?}"),
    }
}

fn list<T>(value: &Value, each: impl Fn(&Value) -> T) -> Vec<T> {
    match value {
        Value::List(items) => items.iter().map(each).collect(),
        other => unreachable!("a SEQUENCE OF of the table, not {other:?}"),
    }
}

fn choice(value: &Value) -> (usize, &Value) {
    match value {
        Value::Choice(index, value) => (usize::from(*index), value),
        other => unreachable!("a CHOICE of the table, not {other:?}"),
    }
}

/// An ENUMERATED value: its name in the list of the member's type, and its number.
fn enumerated(value: &Value, ty: u16, name: &str) -> (&'static str, i32) {
    let Value::Enum(number) = value else {
        unreachable!("an ENUMERATED of the table, not {value:?}")
    };
    let member_ty = members(ty)[member_index(ty, name)].ty;
    let Kind::Enumerated { values, .. } = resolved(member_ty).kind else {
        unreachable!("{name} is an ENUMERATED")
    };
    let name = values
        .iter()
        .find(|(_, value)| value == number)
        .map(|(name, _)| *name)
        .expect("an ENUMERATED number has a name");
    (name, *number)
}

fn enum_number(value: &Value) -> i32 {
    match value {
        Value::Enum(number) => *number,
        other => unreachable!("an ENUMERATED of the table, not {other:?}"),
    }
}

fn object_id(value: &Value) -> ObjectId {
    match choice(value) {
        (0, value) => ObjectId::Id(int4(value)),
        (_, value) => ObjectId::Str(string(value)),
    }
}

/// The model of a `PssmWithParameters` value, or LOSAT's refusal.
fn pssm_of(value: &Value) -> Result<PssmWithParameters, String> {
    const T: u16 = named::PSSMWITHPARAMETERS;
    Ok(PssmWithParameters {
        pssm: pssm(member(value, T, "pssm").expect("mandatory"))?,
        params: member(value, T, "params").map(params),
    })
}

fn pssm(value: &Value) -> Result<Pssm, String> {
    const T: u16 = named::PSSM;
    let get = |name: &str| member(value, T, name);
    Ok(Pssm {
        is_protein: get("isProtein").map(boolean),
        identifier: get("identifier").map(object_id),
        num_rows: get("numRows").map(int4).expect("mandatory"),
        num_columns: get("numColumns").map(int4).expect("mandatory"),
        row_labels: get("rowLabels").map(|labels| list(labels, string)),
        by_row: get("byRow").map(boolean),
        query: get("query").map(seq_entry).transpose()?,
        intermediate_data: get("intermediateData").map(intermediate_data),
        final_data: get("finalData").map(final_data),
    })
}

/// The query's `Seq-entry`: its `seq` variant, a Bioseq.
fn seq_entry(value: &Value) -> Result<PssmQuery, String> {
    let (variant, bioseq) = choice(value);
    let name = members(named::SEQ_ENTRY)[variant].name;
    if name != "seq" {
        return Err(format!(
            "a query Seq-entry {name} is not supported by LOSAT's PSSM reader"
        ));
    }
    const T: u16 = named::BIOSEQ;
    let get = |name: &str| member(bioseq, T, name);
    Ok(PssmQuery {
        ids: list(get("id").expect("mandatory"), |id| match id {
            Value::SeqId(id) => id.clone(),
            other => unreachable!("a Seq-id is read into the Seq-id model, not {other:?}"),
        }),
        descr: get("descr").map(|descr| {
            list(descr, |descriptor| {
                let (variant, value) = choice(descriptor);
                match members(named::SEQDESC)[variant].name {
                    "title" => Seqdesc::Title(string(value)),
                    _ => Seqdesc::Other(descriptor.clone()),
                }
            })
        }),
        inst: seq_inst(get("inst").expect("mandatory"))?,
        annot: get("annot").map(|annot| list(annot, Value::clone)),
    })
}

fn seq_inst(value: &Value) -> Result<SeqInst, String> {
    const T: u16 = named::SEQ_INST;
    let get = |name: &str| member(value, T, name);
    let seq_data = match get("seq-data") {
        None => None,
        Some(data) => {
            let (variant, value) = choice(data);
            Some(match members(named::SEQ_DATA)[variant].name {
                "iupacaa" => SeqData::Iupacaa(string(value)),
                "ncbieaa" => SeqData::Ncbieaa(string(value)),
                "ncbistdaa" => SeqData::Ncbistdaa(match value {
                    Value::Octets(codes) => codes.clone(),
                    other => unreachable!("NCBIstdaa is an OCTET STRING, not {other:?}"),
                }),
                other => {
                    return Err(format!(
                        "query Seq-data {other} is not supported by LOSAT's PSSM reader"
                    ))
                }
            })
        }
    };
    Ok(SeqInst {
        repr: enumerated(get("repr").expect("mandatory"), T, "repr"),
        mol: enumerated(get("mol").expect("mandatory"), T, "mol"),
        length: get("length").map(|length| u32::try_from(int(length)).expect("a Uint4")),
        fuzz: get("fuzz").cloned(),
        topology: get("topology").map(enum_number),
        strand: get("strand").map(enum_number),
        seq_data,
        ext: get("ext").cloned(),
        hist: get("hist").cloned(),
    })
}

fn intermediate_data(value: &Value) -> PssmIntermediateData {
    const T: u16 = named::PSSMINTERMEDIATEDATA;
    let get = |name: &str| member(value, T, name);
    let ints = |name: &str| get(name).map(|values| list(values, int4));
    let reals = |name: &str| get(name).map(|values| list(values, real));
    PssmIntermediateData {
        res_freqs_per_pos: ints("resFreqsPerPos"),
        weighted_res_freqs_per_pos: reals("weightedResFreqsPerPos"),
        freq_ratios: reals("freqRatios").expect("mandatory"),
        information_content: reals("informationContent"),
        gapless_column_weights: reals("gaplessColumnWeights"),
        sigma: reals("sigma"),
        interval_sizes: ints("intervalSizes"),
        num_matching_seqs: ints("numMatchingSeqs"),
        num_indept_obsr: reals("numIndeptObsr"),
    }
}

fn final_data(value: &Value) -> PssmFinalData {
    const T: u16 = named::PSSMFINALDATA;
    let get = |name: &str| member(value, T, name);
    PssmFinalData {
        scores: get("scores")
            .map(|scores| list(scores, int4))
            .expect("mandatory"),
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

fn params(value: &Value) -> PssmParameters {
    const T: u16 = named::PSSMPARAMETERS;
    let get = |name: &str| member(value, T, name);
    PssmParameters {
        pseudocount: get("pseudocount").map(int4),
        rpsdbparams: get("rpsdbparams").map(|rps| {
            const R: u16 = named::FORMATRPSDBPARAMETERS;
            FormatRpsDbParameters {
                matrix_name: member(rps, R, "matrixName").map(string).expect("mandatory"),
                gap_open: member(rps, R, "gapOpen").map(int4),
                gap_extend: member(rps, R, "gapExtend").map(int4),
            }
        }),
        constraints: get("constraints").cloned(),
        bit_score_thresh: get("bitScoreThresh").map(real),
        bit_score_reporting_thresh: get("bitScoreReportingThresh").map(real),
        annotated_sites: get("annotatedSites").cloned(),
    }
}

#[cfg(test)]
#[path = "pssm_text_tests.rs"]
mod tests;
