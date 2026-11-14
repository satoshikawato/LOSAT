//! What NCBI's hand-written classes do when the text reader has read one of them (their
//! `PostRead`), for the classes of the table whose `PostRead` changes the object or stops the
//! reader (SURVEY_D3 §2.4 of line RG series F: `CSeq_descr`, `CVariation_ref`; the others keep
//! caches only), and the Seq-id, which LOSAT keeps in its one Seq-id model (`seqid::model`).

use super::read::{Post, ReadError};
use super::table::named;
use super::{member_index, Value};
use crate::seqid::model::{
    Date, DateStd, Dbtag, GiimportId, IdPat, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId,
    PdbSeqId, TextseqId,
};

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:1497-1527
// ```c++
// void CDense_seg::CReserveHook::PreReadClassMember(CObjectIStream& in,
//                                                   const CObjectInfoMI& member)
// {
//     static CSafeStatic<NCBI_PARAM_TYPE(OBJECTS, DENSE_SEG_RESERVE)> s_Reserve;
//
//     if ( !s_Reserve->Get() ) {
//         return;
//     }
//     CDense_seg& ds = *CType<CDense_seg>::Get(member.GetClassObject());
//     size_t numseg = ds.GetNumseg();
//     try {
//         switch ( member.GetMemberIndex() ) {
//         case 4: // "starts"
//             ds.SetStarts().reserve(ds.GetDim()*numseg);
//             break;
//         case 5: // "lens"
//             ds.SetLens().reserve(numseg);
//             break;
//         case 6: // "strands"
//             ds.SetStrands().reserve(ds.GetDim()*numseg);
//             break;
//         default:
//             break;
//         }
//     }
//     catch ( bad_alloc& /*ignored*/ ) {
//         // ignore insufficient memory exception from advisory reserve()
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/include/objects/seqalign/Dense_seg.hpp:181-182 (the hook is global:
// `NCBISER_HAVE_GLOBAL_READ_MEMBER_HOOK(CDense_seg, "starts,lens,strands", ...)`)
// `numseg` and `dim` are `int` (Dense_seg_.hpp:92-93), so a negative one is a huge `size_t`, and the
// product wraps. `reserve` throws `std::length_error("vector::reserve")` above `max_size()`,
// PTRDIFF_MAX / 4 = 2^61 - 1 for the 4-byte elements of the three vectors (NCBI 2.17.0 probes p5:
// dim 2^30 x numseg 2^31-1 reads on, dim 2^30+1 stops), which no frame catches (they catch
// `CException`s): `CATCH_ALL`'s `std::exception` clause reports it, exit status 255. A smaller
// request that cannot be allocated is a `bad_alloc`, which the hook ignores, as LOSAT does not
// allocate. (`OBJECTS_DENSE_SEG_RESERVE` is NCBI's default, true: LOSAT reads no NCBI parameter.)
// The reserve hooks of Seq-graph (`numval` is a TSeqPos: never above `max_size()`) and of
// SeqTable-multi-data (only in a Seq-annot, which the PSSM search refuses) need nothing here.
/// The hook of NCBI's class `ty` before the value of its member at position `member` is read,
/// with the members `read` so far.
pub(super) fn pre_read_member(
    ty: u16,
    member: usize,
    read: &[(u16, Value)],
) -> Result<(), ReadError> {
    if ty != named::DENSE_SEG {
        return Ok(());
    }
    let get = |name: &str| {
        let index = member_index(named::DENSE_SEG, name);
        read.iter()
            .find(|(position, _)| usize::from(*position) == index)
            .map(|(_, value)| int(value))
    };
    // `GetDim()`: the DEFAULT 2 when the member was not read.
    let dim = get("dim").unwrap_or(2) as i32 as i64 as u64;
    let numseg = get("numseg").unwrap_or(0) as i32 as i64 as u64;
    const MAX_SIZE: u64 = (i64::MAX as u64) / 4;
    let size = match super::members(ty)[member].name {
        "starts" | "strands" => dim.wrapping_mul(numseg),
        "lens" => numseg,
        _ => return Ok(()),
    };
    if size > MAX_SIZE {
        return Err(ReadError::Std("vector::reserve"));
    }
    Ok(())
}

/// The hook of the type `ty` after its value was read (not after a skip).
pub(super) fn post_read(ty: u16, value: Value, posts: &mut Vec<Post>) -> Result<Value, ReadError> {
    match ty {
        named::SEQ_DESCR => seq_descr_post_read(value),
        named::VARIATION_REF => Ok(variation_ref_post_read(value, posts)),
        named::SEQ_ID => Ok(Value::SeqId(seq_id_of(&value))),
        _ => Ok(value),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seq/Seq_descr.cpp:62-68
// ```c++
// void CSeq_descr::PostRead(void) const
// {
//     if ( !s_SeqDescrAllowEmpty->Get() && Get().empty() ) {
//         NCBI_THROW(CSerialException, eInvalidData,
//                    "empty Seq-descr is not allowed");
//     }
// }
// ```
// (`OBJECTS_SEQ_DESCR_ALLOW_EMPTY` is false unless NCBI's environment or configuration file sets
// it, which LOSAT does not read; NCBI 2.17.0 gives the error with the variable set too, c23 of
// SURVEY_D3.)
fn seq_descr_post_read(value: Value) -> Result<Value, ReadError> {
    match &value {
        Value::List(items) if items.is_empty() => {
            Err(ReadError::Serial("empty Seq-descr is not allowed".into()))
        }
        _ => Ok(value),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Variation_ref.cpp:56-151
// ```c++
// void CVariation_ref::PostRead()
// {
//     // population-data: ignore; deprecated, drop
//     if (Tparent::IsSetPopulation_data()) {
//         ERR_POST(Error
//                  << "Variation-ref.population-data is deprecated and "
//                  "will be ignored");
//         Tparent::ResetPopulation_data();
//     }
//
//     // validated: move to VariantProperties.other-validation
//     if (Tparent::IsSetValidated()) {
//         if (SetVariant_prop().IsSetOther_validation()) {
//             ERR_POST(Error
//                      << "Both Variation-ref.validated and "
//                      "Variation-ref.variant-properties.other-validation are "
//                      "set; ignoring Variation-ref.validated");
//         }
//         else {
//             SetVariant_prop().SetOther_validation(Tparent::GetValidated());
//         }
//         Tparent::ResetValidated();
//     }
//
//     // clinical-test: should be moved to Seq-feat.dbxref; no access here!
//     // FIXME: is this used anywhere?
//     if (Tparent::IsSetClinical_test()) {
//         ERR_POST(Error
//                  << "Variation-ref.clinical-test is deprecated and "
//                  "will be ignored");
//         Tparent::ResetClinical_test();
//     }
//
//     // allele-origin: move to Variant-properties
//     if (Tparent::IsSetAllele_origin()) {
//         if (SetVariant_prop().IsSetAllele_origin()) {
//             ERR_POST(Error
//                      << "Both Variation-ref.allele-origin and "
//                      "Variation-ref.variant-properties.allele-origin are "
//                      "set; ignoring Variation-ref.validated");
//         }
//         else {
//             SetVariant_prop().SetAllele_origin(Tparent::GetAllele_origin());
//         }
//         Tparent::ResetAllele_origin();
//     }
//     ... (allele-state, allele-frequency and is-ancestral-allele the same way)
//
//     // pub: move to Seq-feat.dbxref; no access here
//     // FIXME: do we need to do this?
//     if (Tparent::IsSetPub()) {
//         ERR_POST(Error
//                  << "Variation-ref.pub is deprecated and "
//                  "will be ignored");
//         Tparent::ResetPub();
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/VariantProperties.cpp:50-55
// ```c++
// CVariantProperties::CVariantProperties(void)
// {
//     // Per Ming Ward and Lon Phan, 4/21/2011
//     // the current bitfield version is 5
//     SetVersion(5);
// }
// ```
// `SetVariant_prop()` makes a `VariantProperties` (version 5) when the member is not set. The
// messages are posted when the class has been read, before the members that follow it.
fn variation_ref_post_read(mut value: Value, posts: &mut Vec<Post>) -> Value {
    let member = |name: &str| member_index(named::VARIATION_REF, name);
    let drop = |value: &mut Value, posts: &mut Vec<Post>, name: &str, message: &'static str| {
        if value.take_member(member(name)).is_some() {
            posts.push(Post::Error(message));
        }
    };
    let moved = |value: &mut Value,
                 posts: &mut Vec<Post>,
                 name: &str,
                 property: &str,
                 message: &'static str| {
        if let Some(moved) = value.take_member(member(name)) {
            let prop = variant_prop(value);
            let index = member_index(named::VARIANTPROPERTIES, property);
            if prop.member(index).is_some() {
                posts.push(Post::Error(message));
            } else {
                prop.set_member(index, moved);
            }
        }
    };
    drop(
        &mut value,
        posts,
        "population-data",
        "Variation-ref.population-data is deprecated and will be ignored",
    );
    moved(
        &mut value,
        posts,
        "validated",
        "other-validation",
        "Both Variation-ref.validated and Variation-ref.variant-properties.other-validation are set; ignoring Variation-ref.validated",
    );
    drop(
        &mut value,
        posts,
        "clinical-test",
        "Variation-ref.clinical-test is deprecated and will be ignored",
    );
    moved(
        &mut value,
        posts,
        "allele-origin",
        "allele-origin",
        "Both Variation-ref.allele-origin and Variation-ref.variant-properties.allele-origin are set; ignoring Variation-ref.validated",
    );
    moved(
        &mut value,
        posts,
        "allele-state",
        "allele-state",
        "Both Variation-ref.allele-state and Variation-ref.variant-properties.allele-state are set; ignoring Variation-ref.validated",
    );
    moved(
        &mut value,
        posts,
        "allele-frequency",
        "allele-frequency",
        "Both Variation-ref.allele-frequency and Variation-ref.variant-properties.allele-frequency are set; ignoring Variation-ref.validated",
    );
    moved(
        &mut value,
        posts,
        "is-ancestral-allele",
        "is-ancestral-allele",
        "Both Variation-ref.is-ancestral-allele and Variation-ref.variant-properties.is-ancestral-allele are set; ignoring Variation-ref.validated",
    );
    drop(
        &mut value,
        posts,
        "pub",
        "Variation-ref.pub is deprecated and will be ignored",
    );
    value
}

/// `SetVariant_prop()`: the member, made (version 5) when it is not set.
fn variant_prop(value: &mut Value) -> &mut Value {
    let index = member_index(named::VARIATION_REF, "variant-prop");
    if value.member(index).is_none() {
        let version = member_index(named::VARIANTPROPERTIES, "version");
        value.set_member(
            index,
            Value::Sequence(vec![(version as u16, Value::Int(5))]),
        );
    }
    match value {
        Value::Sequence(members) => {
            &mut members
                .iter_mut()
                .find(|(position, _)| usize::from(*position) == index)
                .expect("the member was set")
                .1
        }
        _ => unreachable!("a Variation-ref is a SEQUENCE"),
    }
}

fn int(value: &Value) -> i64 {
    match value {
        Value::Int(n) => *n,
        other => unreachable!("an INTEGER of the table, not {other:?}"),
    }
}

fn int4(value: &Value) -> i32 {
    i32::try_from(int(value)).expect("an Int4 of the table")
}

fn string(value: &Value) -> Vec<u8> {
    match value {
        Value::Str(text) => text.clone(),
        other => unreachable!("a string of the table, not {other:?}"),
    }
}

fn choice(value: &Value) -> (usize, &Value) {
    match value {
        Value::Choice(index, value) => (usize::from(*index), value),
        other => unreachable!("a CHOICE of the table, not {other:?}"),
    }
}

/// The member `name` of the SEQUENCE `value` of the type `ty`, when set.
fn field<'v>(value: &'v Value, ty: u16, name: &str) -> Option<&'v Value> {
    value.member(member_index(ty, name))
}

/// The Seq-id `value` (read under the table) in LOSAT's Seq-id model, which holds each member of
/// every alternative (an OPTIONAL member as `None` when not set, `PDB-seq-id.chain` set to its
/// DEFAULT as `Some(32)`), so the Seq-id keeps what was read.
pub(super) fn seq_id_of(value: &Value) -> ParsedSeqId {
    let (variant, value) = choice(value);
    let textseq = |value: &Value| TextseqId {
        name: field(value, named::TEXTSEQ_ID, "name").map(string),
        accession: field(value, named::TEXTSEQ_ID, "accession").map(string),
        release: field(value, named::TEXTSEQ_ID, "release").map(string),
        version: field(value, named::TEXTSEQ_ID, "version").map(int4),
    };
    match crate::seqid::model::SEQ_ID_ALTERNATIVES[variant] {
        "local" => ParsedSeqId::Local(object_id(value)),
        "gibbsq" => ParsedSeqId::Gibbsq(int4(value)),
        "gibbmt" => ParsedSeqId::Gibbmt(int4(value)),
        "giim" => ParsedSeqId::Giim(GiimportId {
            id: field(value, named::GIIMPORT_ID, "id")
                .map(int4)
                .expect("mandatory"),
            db: field(value, named::GIIMPORT_ID, "db").map(string),
            release: field(value, named::GIIMPORT_ID, "release").map(string),
        }),
        "genbank" => ParsedSeqId::Genbank(textseq(value)),
        "embl" => ParsedSeqId::Embl(textseq(value)),
        "pir" => ParsedSeqId::Pir(textseq(value)),
        "swissprot" => ParsedSeqId::Swissprot(textseq(value)),
        "patent" => {
            let cit = field(value, named::PATENT_SEQ_ID, "cit").expect("mandatory");
            let (number, id) = choice(field(cit, named::ID_PAT, "id").expect("mandatory"));
            ParsedSeqId::Patent(PatentSeqId {
                seqid: field(value, named::PATENT_SEQ_ID, "seqid")
                    .map(int4)
                    .expect("mandatory"),
                cit: IdPat {
                    country: field(cit, named::ID_PAT, "country")
                        .map(string)
                        .expect("mandatory"),
                    id: if number == 0 {
                        PatentNumber::Number(string(id))
                    } else {
                        PatentNumber::AppNumber(string(id))
                    },
                    doc_type: field(cit, named::ID_PAT, "doc-type").map(string),
                },
            })
        }
        "other" => ParsedSeqId::Other(textseq(value)),
        "general" => ParsedSeqId::General(dbtag(value)),
        "gi" => ParsedSeqId::Gi(int(value)),
        "ddbj" => ParsedSeqId::Ddbj(textseq(value)),
        "prf" => ParsedSeqId::Prf(textseq(value)),
        "pdb" => ParsedSeqId::Pdb(PdbSeqId {
            mol: field(value, named::PDB_SEQ_ID, "mol")
                .map(string)
                .expect("mandatory"),
            chain: field(value, named::PDB_SEQ_ID, "chain").map(int4),
            rel: field(value, named::PDB_SEQ_ID, "rel").map(date),
            chain_id: field(value, named::PDB_SEQ_ID, "chain-id").map(string),
        }),
        "tpg" => ParsedSeqId::Tpg(textseq(value)),
        "tpe" => ParsedSeqId::Tpe(textseq(value)),
        "tpd" => ParsedSeqId::Tpd(textseq(value)),
        "gpipe" => ParsedSeqId::Gpipe(textseq(value)),
        _ => ParsedSeqId::NamedAnnotTrack(textseq(value)),
    }
}

fn object_id(value: &Value) -> ObjectId {
    match choice(value) {
        (0, value) => ObjectId::Id(int4(value)),
        (_, value) => ObjectId::Str(string(value)),
    }
}

fn dbtag(value: &Value) -> Dbtag {
    Dbtag {
        db: field(value, named::DBTAG, "db")
            .map(string)
            .expect("mandatory"),
        tag: field(value, named::DBTAG, "tag")
            .map(object_id)
            .expect("mandatory"),
    }
}

fn date(value: &Value) -> Date {
    match choice(value) {
        (0, value) => Date::Str(string(value)),
        (_, value) => {
            let get = |name: &str| field(value, named::DATE_STD, name);
            Date::Std(DateStd {
                year: get("year").map(int4).expect("mandatory"),
                month: get("month").map(int4),
                day: get("day").map(int4),
                season: get("season").map(string),
                hour: get("hour").map(int4),
                minute: get("minute").map(int4),
                second: get("second").map(int4),
            })
        }
    }
}
