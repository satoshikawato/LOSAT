//! The reader of a PSSM in ASN.1 text, ASN.1 binary or XML (`PssmWithParameters`, the file of
//! psiblast's `-out_pssm` and of `-in_pssm`): NCBI's object streams (`CObjectIStreamAsn`,
//! `CObjectIStreamAsnBinary`, `CObjectIStreamXml`) for the types of `report::serial::pssm`, whose
//! writer gives the same bytes back for a PSSM in ASN.1 text that psiblast wrote.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1338-1345
//!
//! The whole document is read by the generic reader of NCBI's type information
//! (`report::serial::asn_types`, line RG item D3-1): every member of every type it reaches, with
//! NCBI's grammar, ranges, replacements and messages (`Post`) and NCBI's errors (`ReadError`). The
//! typed model is then made of the value (`pssm_of`): the members LOSAT's search reads, typed, and
//! the others kept as the generic reader's values: also a query `Seq-entry` that is a `set` and
//! `Seq-data` other than `iupacaa`, `ncbieaa`, `ncbistdaa` and the nucleotide `iupacna`, `ncbi2na`
//! and `ncbi4na` (which makeprofiledb's database writer converts), which NCBI's object manager
//! takes (or stops at) before NCBI reads them; the PSSM search refuses them, and the nucleotide
//! codings, after the object manager's step, for a run that goes on
//! (`pssm_search::in_pssm::psi_search_and_object`). Only the readers' own limits are refused here
//! (`PssmReadError::Refused`).

use crate::report::serial::asn_input::Names;
use crate::report::serial::asn_text_reader::SliceSource;
use crate::report::serial::asn_types::read::ObjectStream;
use crate::report::serial::asn_types::read_xml::XmlStyle;
use crate::report::serial::asn_types::table::{member, named};
use crate::report::serial::asn_types::{
    members, read_ber_document, read_document, resolved, Kind, Post, ReadError, Value, M,
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
    pssm_read(read_document(text, named::PSSMWITHPARAMETERS, Names::PSSM))
}

/// [`read_pssm`] of the `PssmWithParameters` object in ASN.1 binary at the start of `bytes`.
pub(crate) fn read_pssm_binary(bytes: &[u8]) -> PssmRead {
    pssm_read(read_ber_document(
        bytes,
        named::PSSMWITHPARAMETERS,
        Names::PSSM,
    ))
}

/// A `PssmWithParameters` file of makeprofiledb (ASN.1 text, or ASN.1 binary with `-binary`):
/// the same table-driven reader with `Names::PROFILE` (`own_serial_message` false, because
/// makeprofiledb prints NCBI's exception text directly).
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1381-1395
pub(crate) fn read_profile(bytes: &[u8], binary: bool) -> PssmRead {
    pssm_read(if binary {
        read_ber_document(bytes, named::PSSMWITHPARAMETERS, Names::PROFILE)
    } else {
        read_document(bytes, named::PSSMWITHPARAMETERS, Names::PROFILE)
    })
}

/// [`read_pssm`] of the `PssmWithParameters` object in NCBI's XML at the start of `bytes`
/// (the stream read as `MSerial_Xml`, in the ASN.1 style).
pub(crate) fn read_pssm_xml(bytes: &[u8]) -> PssmRead {
    pssm_read(
        ObjectStream::xml(SliceSource::new(bytes), XmlStyle::Asn, Names::PSSM)
            .read_object(named::PSSMWITHPARAMETERS),
    )
}

fn pssm_read((value, posts): (Result<Value, ReadError>, Vec<Post>)) -> PssmRead {
    PssmRead {
        result: value
            .map_err(|error| match error {
                // LOSAT's own limit of a reader: a refusal, not NCBI's error.
                ReadError::Unsupported(what) => PssmReadError::Refused(what),
                error => PssmReadError::Ncbi(error),
            })
            .and_then(pssm_of),
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

/// The member `m` of the SEQUENCE `value`, taken out, when set.
fn take(value: &mut Value, m: M) -> Option<Value> {
    value.take_member(m.index())
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
fn enumerated(value: Value, m: M) -> (&'static str, i32) {
    let number = enum_number(value);
    let Kind::Enumerated { values, .. } = resolved(m.ty).kind else {
        unreachable!("{m:?} is an ENUMERATED")
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
fn pssm_of(mut value: Value) -> Result<PssmWithParameters, PssmReadError> {
    use member::pssmwithparameters as t;
    Ok(PssmWithParameters {
        pssm: pssm(take(&mut value, t::PSSM).expect("mandatory"))?,
        params: take(&mut value, t::PARAMS).map(params),
    })
}

fn pssm(mut value: Value) -> Result<Pssm, PssmReadError> {
    use member::pssm as t;
    let mut get = |m: M| take(&mut value, m);
    Ok(Pssm {
        is_protein: get(t::ISPROTEIN).map(boolean),
        identifier: get(t::IDENTIFIER).map(object_id),
        num_rows: get(t::NUMROWS).map(int4).expect("mandatory"),
        num_columns: get(t::NUMCOLUMNS).map(int4).expect("mandatory"),
        row_labels: get(t::ROWLABELS).map(|labels| list(labels).into_iter().map(string).collect()),
        by_row: get(t::BYROW).map(boolean),
        query: get(t::QUERY).map(seq_entry).transpose()?,
        intermediate_data: get(t::INTERMEDIATEDATA).map(intermediate_data),
        final_data: get(t::FINALDATA).map(final_data),
    })
}

// NCBI reference (598d8ae6): c++/src/objmgr/bioseq_info.cpp:292-308
/// The query's `Seq-entry`: its `seq` variant, a Bioseq, or its `set` variant as read (whose
/// elements stated as `NULL` the object manager's step meets, `objmgr::annot_index`). An
/// element of the Bioseq's `id` stated as `NULL` is refused: NCBI's object manager dereferences
/// each Seq-id when the query joins its scope (`CBioseq_Info`), and its report of the exception
/// prints memory addresses (psiblast and tblastn 2.17.0, rc 255, probe D3c-3b `t_id_elem_null`,
/// `tb_id_elem_null`, `b_id_elem_null`); LOSAT's model holds the Seq-ids typed.
fn seq_entry(value: Value) -> Result<SeqEntry, PssmReadError> {
    let (variant, mut bioseq) = choice(value);
    if members(named::SEQ_ENTRY)[variant].name() != "seq" {
        return Ok(SeqEntry::Set(bioseq));
    }
    use member::bioseq as t;
    let mut get = |m: M| take(&mut bioseq, m);
    let ids = list(get(t::ID).expect("mandatory"))
        .into_iter()
        .map(|id| match id {
            Value::Null => Err(PssmReadError::Refused(NULL_SEQ_ID.to_string())),
            id => Ok(id
                .into_seq_id()
                .unwrap_or_else(|| unreachable!("a Seq-id is read into the Seq-id model"))),
        })
        .collect::<Result<_, _>>()?;
    Ok(SeqEntry::Seq(Box::new(PssmQuery {
        ids,
        descr: get(t::DESCR).map(|descr| {
            list(descr)
                .into_iter()
                .map(|mut descriptor| {
                    match &mut descriptor {
                        Value::Choice(variant, title)
                            if members(named::SEQDESC)[usize::from(*variant)].name() == "title" =>
                        {
                            return Seqdesc::Title(string(title.take()));
                        }
                        // An element stated as `NULL` (`asn_types::Value::Null`).
                        Value::Null => return Seqdesc::Empty,
                        _ => {}
                    }
                    Seqdesc::Other(descriptor)
                })
                .collect()
        }),
        inst: seq_inst(get(t::INST).expect("mandatory")),
        annot: get(t::ANNOT).map(list),
    })))
}

/// The refusal of an element of the query's `id` stated as `NULL` (`seq_entry`).
const NULL_SEQ_ID: &str = "an element of PssmWithParameters.pssm.query.seq.id stated as NULL, \
     which NCBI's object manager dereferences when the query joins its scope (a report of memory \
     addresses)";

fn seq_inst(mut value: Value) -> SeqInst {
    use member::seq_inst as t;
    let mut get = |m: M| take(&mut value, m);
    let seq_data = match get(t::SEQ_DATA) {
        None => None,
        Some(data) => {
            let (variant, mut value) = choice(data);
            let index = u16::try_from(variant).expect("a CHOICE variant");
            let mut octets = |name: &str| match &mut value {
                Value::Octets(codes) => std::mem::take(codes),
                other => unreachable!("{name} is an OCTET STRING, not {other:?}"),
            };
            Some(match members(named::SEQ_DATA)[variant].name() {
                "iupacna" => SeqData::Iupacna(string(value)),
                "iupacaa" => SeqData::Iupacaa(string(value)),
                "ncbi2na" => SeqData::Ncbi2na(octets("NCBI2na")),
                "ncbi4na" => SeqData::Ncbi4na(octets("NCBI4na")),
                "ncbieaa" => SeqData::Ncbieaa(string(value)),
                "ncbistdaa" => SeqData::Ncbistdaa(octets("NCBIstdaa")),
                _ => SeqData::Other(Value::Choice(index, Box::new(value))),
            })
        }
    };
    SeqInst {
        repr: enumerated(get(t::REPR).expect("mandatory"), t::REPR),
        mol: enumerated(get(t::MOL).expect("mandatory"), t::MOL),
        length: get(t::LENGTH).map(|length| u32::try_from(int(length)).expect("a Uint4")),
        fuzz: get(t::FUZZ),
        topology: get(t::TOPOLOGY).map(enum_number),
        strand: get(t::STRAND).map(enum_number),
        seq_data,
        ext: get(t::EXT),
        hist: get(t::HIST),
    }
}

fn intermediate_data(mut value: Value) -> PssmIntermediateData {
    use member::pssmintermediatedata as t;
    let mut get = |m: M| take(&mut value, m);
    PssmIntermediateData {
        res_freqs_per_pos: get(t::RESFREQSPERPOS).map(int4s),
        weighted_res_freqs_per_pos: get(t::WEIGHTEDRESFREQSPERPOS).map(reals),
        freq_ratios: get(t::FREQRATIOS).map(reals).expect("mandatory"),
        information_content: get(t::INFORMATIONCONTENT).map(reals),
        gapless_column_weights: get(t::GAPLESSCOLUMNWEIGHTS).map(reals),
        sigma: get(t::SIGMA).map(reals),
        interval_sizes: get(t::INTERVALSIZES).map(int4s),
        num_matching_seqs: get(t::NUMMATCHINGSEQS).map(int4s),
        num_indept_obsr: get(t::NUMINDEPTOBSR).map(reals),
    }
}

fn final_data(mut value: Value) -> PssmFinalData {
    use member::pssmfinaldata as t;
    let mut get = |m: M| take(&mut value, m);
    PssmFinalData {
        scores: get(t::SCORES).map(int4s).expect("mandatory"),
        lambda: get(t::LAMBDA).map(real).expect("mandatory"),
        kappa: get(t::KAPPA).map(real).expect("mandatory"),
        h: get(t::H).map(real).expect("mandatory"),
        scaling_factor: get(t::SCALINGFACTOR).map(int4),
        lambda_ungapped: get(t::LAMBDAUNGAPPED).map(real),
        kappa_ungapped: get(t::KAPPAUNGAPPED).map(real),
        h_ungapped: get(t::HUNGAPPED).map(real),
        word_score_threshold: get(t::WORDSCORETHRESHOLD).map(real),
    }
}

fn params(mut value: Value) -> PssmParameters {
    use member::pssmparameters as t;
    let mut get = |m: M| take(&mut value, m);
    PssmParameters {
        pseudocount: get(t::PSEUDOCOUNT).map(int4),
        rpsdbparams: get(t::RPSDBPARAMS).map(|mut rps| {
            use member::formatrpsdbparameters as r;
            FormatRpsDbParameters {
                matrix_name: take(&mut rps, r::MATRIXNAME)
                    .map(string)
                    .expect("mandatory"),
                gap_open: take(&mut rps, r::GAPOPEN).map(int4),
                gap_extend: take(&mut rps, r::GAPEXTEND).map(int4),
            }
        }),
        constraints: get(t::CONSTRAINTS),
        bit_score_thresh: get(t::BITSCORETHRESH).map(real),
        bit_score_reporting_thresh: get(t::BITSCOREREPORTINGTHRESH).map(real),
        annotated_sites: get(t::ANNOTATEDSITES),
    }
}

#[cfg(test)]
#[path = "pssm_text_tests.rs"]
mod tests;
