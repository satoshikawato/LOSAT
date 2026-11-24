//! The reader of the Seq-id of NCBI's ASN.1 text (`CObjectIStreamAsn` of a `CSeq_id`), every
//! alternative, as `ber` is its BER reader: the rpsbproc archive reader reads its Seq-ids with it.
//! Each value is read below the frame `at` that holds it (`asn_text_reader::Frame`), and the
//! posts of its strings are given to `posts` as they are read (`Reader::string`).
//! (The PSSM reader, psiblast and tblastn `-in_pssm`, reads them with the generic typed reader,
//! `report::serial::asn_types`, and converts them to the same model.)

use crate::report::serial::asn_text_reader::{AsnSource, Frame, Reader, ReplacedSink, Result};
use crate::seqid::model::{
    Date, DateStd, Dbtag, GiimportId, IdPat, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId,
    PdbSeqId, TextseqId, SEQ_ID_ALTERNATIVES,
};

impl<S: AsnSource> Reader<S> {
    // NCBI reference (598d8ae6): c++/src/objects/general/general.asn:78-80
    pub(crate) fn object_id<K: ReplacedSink + ?Sized>(
        &mut self,
        at: &Frame<'_>,
        posts: &mut K,
    ) -> Result<ObjectId> {
        let variant = self.choice(at, &["id", "str"])?;
        Ok(if variant.index() == 0 {
            ObjectId::Id(self.int4()?)
        } else {
            ObjectId::Str(self.string(&variant, posts)?)
        })
    }

    // NCBI reference (598d8ae6): c++/src/objects/general/general.asn:71-73
    pub(crate) fn dbtag<K: ReplacedSink + ?Sized>(
        &mut self,
        at: &Frame<'_>,
        posts: &mut K,
    ) -> Result<Dbtag> {
        let mut members = self.begin_members(at, &["db", "tag"])?;
        let (mut db, mut tag) = (None, None);
        while let Some(index) = self.next_member(&mut members)? {
            match index {
                0 => db = Some(self.string(members.at(), posts)?),
                _ => tag = Some(self.object_id(members.at(), posts)?),
            }
        }
        match (db, tag) {
            (Some(db), Some(tag)) => Ok(Dbtag { db, tag }),
            _ => Err(self.error("a mandatory member is missing")),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objects/general/general.asn:54-65
    pub(crate) fn date<K: ReplacedSink + ?Sized>(
        &mut self,
        at: &Frame<'_>,
        posts: &mut K,
    ) -> Result<Date> {
        let variant = self.choice(at, &["str", "std"])?;
        if variant.index() == 0 {
            return Ok(Date::Str(self.string(&variant, posts)?));
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
        let mut members = self.begin_members(&variant, &NAMES)?;
        while let Some(index) = self.next_member(&mut members)? {
            match index {
                0 => {
                    date.year = self.int4()?;
                    year = true;
                }
                1 => date.month = Some(self.int4()?),
                2 => date.day = Some(self.int4()?),
                3 => date.season = Some(self.string(members.at(), posts)?),
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
    pub(crate) fn seq_id<K: ReplacedSink + ?Sized>(
        &mut self,
        at: &Frame<'_>,
        posts: &mut K,
    ) -> Result<ParsedSeqId> {
        let at = &self.choice(at, &SEQ_ID_ALTERNATIVES)?;
        let variant = at.index();
        Ok(match variant {
            0 => ParsedSeqId::Local(self.object_id(at, posts)?),
            1 => ParsedSeqId::Gibbsq(self.int4()?),
            2 => ParsedSeqId::Gibbmt(self.int4()?),
            3 => {
                let mut members = self.begin_members(at, &["id", "db", "release"])?;
                let mut giim = GiimportId {
                    id: 0,
                    db: None,
                    release: None,
                };
                let mut id = false;
                while let Some(index) = self.next_member(&mut members)? {
                    match index {
                        0 => {
                            giim.id = self.int4()?;
                            id = true;
                        }
                        1 => giim.db = Some(self.string(members.at(), posts)?),
                        _ => giim.release = Some(self.string(members.at(), posts)?),
                    }
                }
                if !id {
                    return Err(self.error("member id missing"));
                }
                ParsedSeqId::Giim(giim)
            }
            8 => {
                let mut members = self.begin_members(at, &["seqid", "cit"])?;
                let (mut seqid, mut cit) = (None, None);
                while let Some(index) = self.next_member(&mut members)? {
                    if index == 0 {
                        seqid = Some(self.int4()?);
                        continue;
                    }
                    let mut members =
                        self.begin_members(members.at(), &["country", "id", "doc-type"])?;
                    let (mut country, mut number, mut doc_type) = (None, None, None);
                    while let Some(index) = self.next_member(&mut members)? {
                        match index {
                            0 => country = Some(self.string(members.at(), posts)?),
                            1 => {
                                let variant =
                                    self.choice(members.at(), &["number", "app-number"])?;
                                let app = variant.index() == 1;
                                let text = self.string(&variant, posts)?;
                                number = Some(if app {
                                    PatentNumber::AppNumber(text)
                                } else {
                                    PatentNumber::Number(text)
                                });
                            }
                            _ => doc_type = Some(self.string(members.at(), posts)?),
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
            10 => ParsedSeqId::General(self.dbtag(at, posts)?),
            11 => ParsedSeqId::Gi(self.int8()?),
            14 => {
                let mut members = self.begin_members(at, &["mol", "chain", "rel", "chain-id"])?;
                let mut pdb = PdbSeqId {
                    mol: Vec::new(),
                    chain: None,
                    rel: None,
                    chain_id: None,
                };
                let mut mol = false;
                while let Some(index) = self.next_member(&mut members)? {
                    match index {
                        0 => {
                            pdb.mol = self.string(members.at(), posts)?;
                            mol = true;
                        }
                        1 => pdb.chain = Some(self.int4()?),
                        2 => pdb.rel = Some(self.date(members.at(), posts)?),
                        _ => pdb.chain_id = Some(self.string(members.at(), posts)?),
                    }
                }
                if !mol {
                    return Err(self.error("member mol missing"));
                }
                ParsedSeqId::Pdb(pdb)
            }
            _ => {
                let mut members =
                    self.begin_members(at, &["name", "accession", "release", "version"])?;
                let mut id = TextseqId::default();
                while let Some(index) = self.next_member(&mut members)? {
                    match index {
                        0 => id.name = Some(self.string(members.at(), posts)?),
                        1 => id.accession = Some(self.string(members.at(), posts)?),
                        2 => id.release = Some(self.string(members.at(), posts)?),
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
