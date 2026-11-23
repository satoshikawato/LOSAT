//! The Seq-id read into LOSAT's Seq-id model (`seqid::model::ParsedSeqId`, which holds each
//! member of every alternative: an OPTIONAL or DEFAULT member as `None` when the stream does not
//! state it, as NCBI keeps its set flag), on the object stack of `read` (`super::typed`): the one
//! reader of an ASN.1 Seq-id (text and binary) of the table's types, for the generic value
//! (`Value::SeqId`) and for makeblastdb's mask files (`writedb::mask_info`).

use super::read::TypedReader;
use super::table::named;
use crate::report::serial::asn_input::{Input, R};
use crate::seqid::model::{
    Date, DateStd, Dbtag, GiimportId, IdPat, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId,
    PdbSeqId, TextseqId,
};

/// A `Seq-id` of any alternative (the variants in the order of `seqid::model::SEQ_ID_ALTERNATIVES`,
/// checked by `seq_id_tests`).
pub(crate) fn read_seq_id<I: Input>(r: &mut TypedReader<I>) -> R<ParsedSeqId> {
    let mut seq_id = None;
    read_seq_id_into(r, &mut seq_id)?;
    Ok(mandatory(seq_id))
}

/// `read_seq_id` into `slot`, whose earlier value it replaces: a reader of many Seq-ids (the
/// masks of a mask file) keeps one slot and moves no Seq-id through the frames.
pub(crate) fn read_seq_id_into<I: Input>(
    r: &mut TypedReader<I>,
    slot: &mut Option<ParsedSeqId>,
) -> R<()> {
    r.choice(
        named::SEQ_ID,
        #[inline(always)]
        |r, index, ty| {
            *slot = Some(match index {
                0 => ParsedSeqId::Local(object_id(r, ty)?),
                1 => ParsedSeqId::Gibbsq(int4(r, ty)?),
                2 => ParsedSeqId::Gibbmt(int4(r, ty)?),
                3 => ParsedSeqId::Giim(giimport_id(r, ty)?),
                4 => ParsedSeqId::Genbank(textseq_id(r, ty)?),
                5 => ParsedSeqId::Embl(textseq_id(r, ty)?),
                6 => ParsedSeqId::Pir(textseq_id(r, ty)?),
                7 => ParsedSeqId::Swissprot(textseq_id(r, ty)?),
                8 => ParsedSeqId::Patent(patent_seq_id(r, ty)?),
                9 => ParsedSeqId::Other(textseq_id(r, ty)?),
                10 => ParsedSeqId::General(dbtag(r, ty)?),
                11 => ParsedSeqId::Gi(r.integer(ty)?),
                12 => ParsedSeqId::Ddbj(textseq_id(r, ty)?),
                13 => ParsedSeqId::Prf(textseq_id(r, ty)?),
                14 => ParsedSeqId::Pdb(pdb_seq_id(r, ty)?),
                15 => ParsedSeqId::Tpg(textseq_id(r, ty)?),
                16 => ParsedSeqId::Tpe(textseq_id(r, ty)?),
                17 => ParsedSeqId::Tpd(textseq_id(r, ty)?),
                18 => ParsedSeqId::Gpipe(textseq_id(r, ty)?),
                _ => ParsedSeqId::NamedAnnotTrack(textseq_id(r, ty)?),
            });
            Ok(())
        },
    )
}

/// An INTEGER of Int4 storage (`ReadInt4`: the number fits).
#[inline]
fn int4<I: Input>(r: &mut TypedReader<I>, ty: u16) -> R<i32> {
    Ok(r.integer(ty)? as i32)
}

/// The value of a member that `ReadClassSequential` read: a member that is neither OPTIONAL nor
/// has a DEFAULT was read when `TypedReader::class` returns.
pub(crate) fn mandatory<T>(value: Option<T>) -> T {
    value.expect("ReadClassSequential has read every mandatory member")
}

/// An `Object-id`.
pub(crate) fn object_id<I: Input>(r: &mut TypedReader<I>, ty: u16) -> R<ObjectId> {
    r.choice(
        ty,
        #[inline(always)]
        |r, index, ty| {
            Ok(if index == 0 {
                ObjectId::Id(int4(r, ty)?)
            } else {
                ObjectId::Str(r.text(ty)?)
            })
        },
    )
}

fn dbtag<I: Input>(r: &mut TypedReader<I>, ty: u16) -> R<Dbtag> {
    let (mut db, mut tag) = (None, None);
    r.class(ty, |r, index, ty| {
        match index {
            0 => db = Some(r.text(ty)?),
            _ => tag = Some(object_id(r, ty)?),
        }
        Ok(())
    })?;
    Ok(Dbtag {
        db: mandatory(db),
        tag: mandatory(tag),
    })
}

fn textseq_id<I: Input>(r: &mut TypedReader<I>, ty: u16) -> R<TextseqId> {
    let mut id = TextseqId::default();
    r.class(ty, |r, index, ty| {
        match index {
            0 => id.name = Some(r.text(ty)?),
            1 => id.accession = Some(r.text(ty)?),
            2 => id.release = Some(r.text(ty)?),
            _ => id.version = Some(int4(r, ty)?),
        }
        Ok(())
    })?;
    Ok(id)
}

fn giimport_id<I: Input>(r: &mut TypedReader<I>, ty: u16) -> R<GiimportId> {
    let (mut id, mut db, mut release) = (None, None, None);
    r.class(ty, |r, index, ty| {
        match index {
            0 => id = Some(int4(r, ty)?),
            1 => db = Some(r.text(ty)?),
            _ => release = Some(r.text(ty)?),
        }
        Ok(())
    })?;
    Ok(GiimportId {
        id: mandatory(id),
        db,
        release,
    })
}

fn patent_seq_id<I: Input>(r: &mut TypedReader<I>, ty: u16) -> R<PatentSeqId> {
    let (mut seqid, mut cit) = (None, None);
    r.class(ty, |r, index, ty| {
        match index {
            0 => seqid = Some(int4(r, ty)?),
            _ => cit = Some(id_pat(r, ty)?),
        }
        Ok(())
    })?;
    Ok(PatentSeqId {
        seqid: mandatory(seqid),
        cit: mandatory(cit),
    })
}

fn id_pat<I: Input>(r: &mut TypedReader<I>, ty: u16) -> R<IdPat> {
    let (mut country, mut id, mut doc_type) = (None, None, None);
    r.class(ty, |r, index, ty| {
        match index {
            0 => country = Some(r.text(ty)?),
            1 => {
                id = Some(r.choice(ty, |r, variant, ty| {
                    Ok(if variant == 0 {
                        PatentNumber::Number(r.text(ty)?)
                    } else {
                        PatentNumber::AppNumber(r.text(ty)?)
                    })
                })?);
            }
            _ => doc_type = Some(r.text(ty)?),
        }
        Ok(())
    })?;
    Ok(IdPat {
        country: mandatory(country),
        id: mandatory(id),
        doc_type,
    })
}

fn pdb_seq_id<I: Input>(r: &mut TypedReader<I>, ty: u16) -> R<PdbSeqId> {
    let (mut mol, mut chain, mut rel, mut chain_id) = (None, None, None, None);
    r.class(ty, |r, index, ty| {
        match index {
            0 => mol = Some(r.text(ty)?),
            1 => chain = Some(int4(r, ty)?),
            2 => rel = Some(date(r, ty)?),
            _ => chain_id = Some(r.text(ty)?),
        }
        Ok(())
    })?;
    Ok(PdbSeqId {
        mol: mandatory(mol),
        chain,
        rel,
        chain_id,
    })
}

fn date<I: Input>(r: &mut TypedReader<I>, ty: u16) -> R<Date> {
    r.choice(ty, |r, variant, ty| {
        if variant == 0 {
            return Ok(Date::Str(r.text(ty)?));
        }
        let mut std = DateStd {
            year: 0,
            month: None,
            day: None,
            season: None,
            hour: None,
            minute: None,
            second: None,
        };
        let mut year = None;
        r.class(ty, |r, index, ty| {
            match index {
                0 => year = Some(int4(r, ty)?),
                1 => std.month = Some(int4(r, ty)?),
                2 => std.day = Some(int4(r, ty)?),
                3 => std.season = Some(r.text(ty)?),
                4 => std.hour = Some(int4(r, ty)?),
                5 => std.minute = Some(int4(r, ty)?),
                _ => std.second = Some(int4(r, ty)?),
            }
            Ok(())
        })?;
        std.year = mandatory(year);
        Ok(Date::Std(std))
    })
}

/// The member names of the types read here in the order the readers above take them: the
/// table's module order (`seq_id_tests`).
#[cfg(test)]
pub(crate) const MEMBER_ORDER: [(u16, &[&str]); 9] = [
    (named::OBJECT_ID, &["id", "str"]),
    (named::DBTAG, &["db", "tag"]),
    (
        named::TEXTSEQ_ID,
        &["name", "accession", "release", "version"],
    ),
    (named::GIIMPORT_ID, &["id", "db", "release"]),
    (named::PATENT_SEQ_ID, &["seqid", "cit"]),
    (named::ID_PAT, &["country", "id", "doc-type"]),
    (named::PDB_SEQ_ID, &["mol", "chain", "rel", "chain-id"]),
    (named::DATE, &["str", "std"]),
    (
        named::DATE_STD,
        &["year", "month", "day", "season", "hour", "minute", "second"],
    ),
];
