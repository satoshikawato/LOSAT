//! The reader of the Seq-id of NCBI's ASN.1 text (`CObjectIStreamAsn` of a `CSeq_id`), every
//! alternative, as `ber` is its BER reader: the rpsbproc archive reader reads its Seq-ids with it.
//! (The PSSM reader, psiblast and tblastn `-in_pssm`, reads them with the generic typed reader,
//! `report::serial::asn_types`, and converts them to the same model.)

use crate::report::serial::asn_text_reader::{AsnSource, Reader, Result};
use crate::seqid::model::{
    Date, DateStd, Dbtag, GiimportId, IdPat, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId,
    PdbSeqId, TextseqId, SEQ_ID_ALTERNATIVES,
};

impl<S: AsnSource> Reader<S> {
    // NCBI reference (598d8ae6): c++/src/objects/general/general.asn:78-80
    pub(crate) fn object_id(&mut self) -> Result<ObjectId> {
        Ok(if self.choice(&["id", "str"])? == 0 {
            ObjectId::Id(self.int4()?)
        } else {
            ObjectId::Str(self.string()?)
        })
    }

    // NCBI reference (598d8ae6): c++/src/objects/general/general.asn:71-73
    pub(crate) fn dbtag(&mut self) -> Result<Dbtag> {
        let mut members = self.begin_members()?;
        let (mut db, mut tag) = (None, None);
        while let Some(index) = self.next_member(&mut members, &["db", "tag"])? {
            match index {
                0 => db = Some(self.string()?),
                _ => tag = Some(self.object_id()?),
            }
        }
        match (db, tag) {
            (Some(db), Some(tag)) => Ok(Dbtag { db, tag }),
            _ => Err(self.error("a mandatory member is missing")),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objects/general/general.asn:54-65
    pub(crate) fn date(&mut self) -> Result<Date> {
        if self.choice(&["str", "std"])? == 0 {
            return Ok(Date::Str(self.string()?));
        }
        const NAMES: [&str; 7] = ["year", "month", "day", "season", "hour", "minute", "second"];
        let mut date = DateStd {
            year: 0,
            month: None,
            day: None,
            season: None,
            hour: None,
            minute: None,
            second: None,
        };
        let mut year = false;
        let mut members = self.begin_members()?;
        while let Some(index) = self.next_member(&mut members, &NAMES)? {
            match index {
                0 => {
                    date.year = self.int4()?;
                    year = true;
                }
                1 => date.month = Some(self.int4()?),
                2 => date.day = Some(self.int4()?),
                3 => date.season = Some(self.string()?),
                4 => date.hour = Some(self.int4()?),
                5 => date.minute = Some(self.int4()?),
                _ => date.second = Some(self.int4()?),
            }
        }
        if !year {
            return Err(self.error("member year missing"));
        }
        Ok(Date::Std(date))
    }

    /// A `Seq-id` of any variant.
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:24-74
    // NCBI reference (598d8ae6): c++/src/objects/biblio/biblio.asn:119-124
    // `gi` is `TGi`, an `Int8` in NCBI's 64-bit build.
    pub(crate) fn seq_id(&mut self) -> Result<ParsedSeqId> {
        let variant = self.choice(&SEQ_ID_ALTERNATIVES)?;
        Ok(match variant {
            0 => ParsedSeqId::Local(self.object_id()?),
            1 => ParsedSeqId::Gibbsq(self.int4()?),
            2 => ParsedSeqId::Gibbmt(self.int4()?),
            3 => {
                let mut members = self.begin_members()?;
                let mut giim = GiimportId {
                    id: 0,
                    db: None,
                    release: None,
                };
                let mut id = false;
                while let Some(index) = self.next_member(&mut members, &["id", "db", "release"])? {
                    match index {
                        0 => {
                            giim.id = self.int4()?;
                            id = true;
                        }
                        1 => giim.db = Some(self.string()?),
                        _ => giim.release = Some(self.string()?),
                    }
                }
                if !id {
                    return Err(self.error("member id missing"));
                }
                ParsedSeqId::Giim(giim)
            }
            8 => {
                let mut members = self.begin_members()?;
                let (mut seqid, mut cit) = (None, None);
                while let Some(index) = self.next_member(&mut members, &["seqid", "cit"])? {
                    if index == 0 {
                        seqid = Some(self.int4()?);
                        continue;
                    }
                    let mut members = self.begin_members()?;
                    let (mut country, mut number, mut doc_type) = (None, None, None);
                    while let Some(index) =
                        self.next_member(&mut members, &["country", "id", "doc-type"])?
                    {
                        match index {
                            0 => country = Some(self.string()?),
                            1 => {
                                let app = self.choice(&["number", "app-number"])? == 1;
                                let text = self.string()?;
                                number = Some(if app {
                                    PatentNumber::AppNumber(text)
                                } else {
                                    PatentNumber::Number(text)
                                });
                            }
                            _ => doc_type = Some(self.string()?),
                        }
                    }
                    let (Some(country), Some(id)) = (country, number) else {
                        return Err(self.error("a mandatory member is missing"));
                    };
                    cit = Some(IdPat {
                        country,
                        id,
                        doc_type,
                    });
                }
                let (Some(seqid), Some(cit)) = (seqid, cit) else {
                    return Err(self.error("a mandatory member is missing"));
                };
                ParsedSeqId::Patent(PatentSeqId { seqid, cit })
            }
            10 => ParsedSeqId::General(self.dbtag()?),
            11 => ParsedSeqId::Gi(self.int8()?),
            14 => {
                let mut members = self.begin_members()?;
                let mut pdb = PdbSeqId {
                    mol: Vec::new(),
                    chain: None,
                    rel: None,
                    chain_id: None,
                };
                let mut mol = false;
                while let Some(index) =
                    self.next_member(&mut members, &["mol", "chain", "rel", "chain-id"])?
                {
                    match index {
                        0 => {
                            pdb.mol = self.string()?;
                            mol = true;
                        }
                        1 => pdb.chain = Some(self.int4()?),
                        2 => pdb.rel = Some(self.date()?),
                        _ => pdb.chain_id = Some(self.string()?),
                    }
                }
                if !mol {
                    return Err(self.error("member mol missing"));
                }
                ParsedSeqId::Pdb(pdb)
            }
            _ => {
                let mut members = self.begin_members()?;
                let mut id = TextseqId::default();
                while let Some(index) =
                    self.next_member(&mut members, &["name", "accession", "release", "version"])?
                {
                    match index {
                        0 => id.name = Some(self.string()?),
                        1 => id.accession = Some(self.string()?),
                        2 => id.release = Some(self.string()?),
                        _ => id.version = Some(self.int4()?),
                    }
                }
                match variant {
                    4 => ParsedSeqId::Genbank(id),
                    5 => ParsedSeqId::Embl(id),
                    6 => ParsedSeqId::Pir(id),
                    7 => ParsedSeqId::Swissprot(id),
                    9 => ParsedSeqId::Other(id),
                    12 => ParsedSeqId::Ddbj(id),
                    13 => ParsedSeqId::Prf(id),
                    15 => ParsedSeqId::Tpg(id),
                    16 => ParsedSeqId::Tpe(id),
                    17 => ParsedSeqId::Tpd(id),
                    18 => ParsedSeqId::Gpipe(id),
                    _ => ParsedSeqId::NamedAnnotTrack(id),
                }
            }
        })
    }
}
