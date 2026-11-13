//! The mask reader builds a `Seq-id` in its own walk (text and binary share one builder,
//! `reader::read_seq_id`) and `seqid::ber::decode_seq_id` builds one for the defline of a
//! database: two readers of the same type. For a valid binary `Seq-id` of every alternative
//! they must give the same `ParsedSeqId`, so that a mask id and a database id compare as NCBI's
//! one `CSeq_id` does.

use super::reader::seq_id_from_binary;
use crate::seqid::ber::{decode_seq_id, Reader};
use crate::seqid::model::{
    Date, DateStd, Dbtag, GiimportId, IdPat, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId,
    PdbSeqId, TextseqId,
};

/// ASN.1 binary with automatic tags: the constructed items NCBI writes (indefinite length,
/// closed by `00 00`) or the same with definite lengths.
struct Writer {
    indefinite: bool,
}

impl Writer {
    fn length(len: usize) -> Vec<u8> {
        if len < 0x80 {
            return vec![len as u8];
        }
        let bytes: Vec<u8> = len
            .to_be_bytes()
            .into_iter()
            .skip_while(|&b| b == 0)
            .collect();
        let mut out = vec![0x80 | bytes.len() as u8];
        out.extend(bytes);
        out
    }

    /// A primitive item: tag, length, content.
    fn primitive(tag: u8, content: &[u8]) -> Vec<u8> {
        let mut out = vec![tag];
        out.extend(Self::length(content.len()));
        out.extend_from_slice(content);
        out
    }

    /// A constructed item.
    fn constructed(&self, tag: u8, content: &[u8]) -> Vec<u8> {
        let mut out = vec![tag];
        if self.indefinite {
            out.push(0x80);
            out.extend_from_slice(content);
            out.extend([0, 0]);
        } else {
            out.extend(Self::length(content.len()));
            out.extend_from_slice(content);
        }
        out
    }

    fn integer(value: i64) -> Vec<u8> {
        let bytes = value.to_be_bytes();
        let mut start = 0;
        while start < 7
            && ((bytes[start] == 0 && bytes[start + 1] & 0x80 == 0)
                || (bytes[start] == 0xff && bytes[start + 1] & 0x80 != 0))
        {
            start += 1;
        }
        Self::primitive(0x02, &bytes[start..])
    }

    fn string(value: &[u8]) -> Vec<u8> {
        Self::primitive(0x1a, value)
    }

    /// `[n] value`: a member of a SEQUENCE or the variant of a CHOICE.
    fn tagged(&self, n: u8, value: &[u8]) -> Vec<u8> {
        self.constructed(0xa0 | n, value)
    }

    fn sequence(&self, members: &[Vec<u8>]) -> Vec<u8> {
        self.constructed(0x30, &members.concat())
    }

    fn object_id(&self, id: &ObjectId) -> Vec<u8> {
        match id {
            ObjectId::Id(n) => self.tagged(0, &Self::integer(i64::from(*n))),
            ObjectId::Str(s) => self.tagged(1, &Self::string(s)),
        }
    }

    fn textseq(&self, id: &TextseqId) -> Vec<u8> {
        let mut members = Vec::new();
        if let Some(name) = &id.name {
            members.push(self.tagged(0, &Self::string(name)));
        }
        if let Some(accession) = &id.accession {
            members.push(self.tagged(1, &Self::string(accession)));
        }
        if let Some(release) = &id.release {
            members.push(self.tagged(2, &Self::string(release)));
        }
        if let Some(version) = id.version {
            members.push(self.tagged(3, &Self::integer(i64::from(version))));
        }
        self.sequence(&members)
    }

    fn date(&self, date: &Date) -> Vec<u8> {
        match date {
            Date::Str(s) => self.tagged(0, &Self::string(s)),
            Date::Std(d) => {
                let int = |n, v: Option<i32>| v.map(|v| self.tagged(n, &Self::integer(v.into())));
                let mut members = vec![self.tagged(0, &Self::integer(d.year.into()))];
                members.extend(int(1, d.month));
                members.extend(int(2, d.day));
                if let Some(season) = &d.season {
                    members.push(self.tagged(3, &Self::string(season)));
                }
                members.extend(int(4, d.hour));
                members.extend(int(5, d.minute));
                members.extend(int(6, d.second));
                self.tagged(1, &self.sequence(&members))
            }
        }
    }

    fn seq_id(&self, id: &ParsedSeqId) -> Vec<u8> {
        use ParsedSeqId as P;
        let (alternative, value) = match id {
            P::Local(o) => (0, self.object_id(o)),
            P::Gibbsq(n) => (1, Self::integer(i64::from(*n))),
            P::Gibbmt(n) => (2, Self::integer(i64::from(*n))),
            P::Giim(g) => {
                let mut members = vec![self.tagged(0, &Self::integer(g.id.into()))];
                if let Some(db) = &g.db {
                    members.push(self.tagged(1, &Self::string(db)));
                }
                if let Some(release) = &g.release {
                    members.push(self.tagged(2, &Self::string(release)));
                }
                (3, self.sequence(&members))
            }
            P::Genbank(t) => (4, self.textseq(t)),
            P::Embl(t) => (5, self.textseq(t)),
            P::Pir(t) => (6, self.textseq(t)),
            P::Swissprot(t) => (7, self.textseq(t)),
            P::Patent(p) => {
                let number = match &p.cit.id {
                    PatentNumber::Number(s) => self.tagged(0, &Self::string(s)),
                    PatentNumber::AppNumber(s) => self.tagged(1, &Self::string(s)),
                };
                let mut cit = vec![
                    self.tagged(0, &Self::string(&p.cit.country)),
                    self.tagged(1, &number),
                ];
                if let Some(doc_type) = &p.cit.doc_type {
                    cit.push(self.tagged(2, &Self::string(doc_type)));
                }
                let members = [
                    self.tagged(0, &Self::integer(p.seqid.into())),
                    self.tagged(1, &self.sequence(&cit)),
                ];
                (8, self.sequence(&members))
            }
            P::Other(t) => (9, self.textseq(t)),
            P::General(d) => {
                let members = [
                    self.tagged(0, &Self::string(&d.db)),
                    self.tagged(1, &self.object_id(&d.tag)),
                ];
                (10, self.sequence(&members))
            }
            P::Gi(n) => (11, Self::integer(*n)),
            P::Ddbj(t) => (12, self.textseq(t)),
            P::Prf(t) => (13, self.textseq(t)),
            P::Pdb(p) => {
                let mut members = vec![self.tagged(0, &Self::string(&p.mol))];
                if let Some(chain) = p.chain {
                    members.push(self.tagged(1, &Self::integer(chain.into())));
                }
                if let Some(rel) = &p.rel {
                    members.push(self.tagged(2, &self.date(rel)));
                }
                if let Some(chain_id) = &p.chain_id {
                    members.push(self.tagged(3, &Self::string(chain_id)));
                }
                (14, self.sequence(&members))
            }
            P::Tpg(t) => (15, self.textseq(t)),
            P::Tpe(t) => (16, self.textseq(t)),
            P::Tpd(t) => (17, self.textseq(t)),
            P::Gpipe(t) => (18, self.textseq(t)),
            P::NamedAnnotTrack(t) => (19, self.textseq(t)),
        };
        self.tagged(alternative, &value)
    }
}

fn b(s: &str) -> Vec<u8> {
    s.as_bytes().to_vec()
}

fn text(name: Option<&str>, acc: Option<&str>, rel: Option<&str>, ver: Option<i32>) -> TextseqId {
    TextseqId {
        name: name.map(b),
        accession: acc.map(b),
        release: rel.map(b),
        version: ver,
    }
}

/// Valid Seq-ids: every alternative, with all members and with the fewest, and numbers that
/// need one, two, three and eight bytes.
fn seq_ids() -> Vec<ParsedSeqId> {
    use ParsedSeqId as P;
    let full = || text(Some("NAME"), Some("AY123456"), Some("rel"), Some(2));
    let least = || text(None, None, None, None);
    let mut ids = vec![
        P::Local(ObjectId::Id(0)),
        P::Local(ObjectId::Id(127)),
        P::Local(ObjectId::Id(128)),
        P::Local(ObjectId::Id(-1)),
        P::Local(ObjectId::Id(-129)),
        P::Local(ObjectId::Id(i32::MIN)),
        P::Local(ObjectId::Id(i32::MAX)),
        P::Local(ObjectId::Str(b("abc"))),
        P::Local(ObjectId::Str(Vec::new())),
        P::Gibbsq(5),
        P::Gibbmt(-70000),
        P::Giim(GiimportId {
            id: 9,
            db: None,
            release: None,
        }),
        P::Giim(GiimportId {
            id: 9,
            db: Some(b("db")),
            release: Some(b("r1")),
        }),
        P::Patent(PatentSeqId {
            seqid: 3,
            cit: IdPat {
                country: b("US"),
                id: PatentNumber::Number(b("123")),
                doc_type: None,
            },
        }),
        P::Patent(PatentSeqId {
            seqid: 3,
            cit: IdPat {
                country: b("EP"),
                id: PatentNumber::AppNumber(b("456")),
                doc_type: Some(b("A1")),
            },
        }),
        P::General(Dbtag {
            db: b("db"),
            tag: ObjectId::Id(12),
        }),
        P::General(Dbtag {
            db: b("db"),
            tag: ObjectId::Str(b("t 1")),
        }),
        P::Gi(0),
        P::Gi(123),
        P::Gi(40_000),
        P::Gi(70_000),
        P::Gi(i64::from(i32::MAX) + 1),
        P::Gi(i64::MAX),
        P::Gi(i64::MIN),
        P::Pdb(PdbSeqId {
            mol: b("1ABC"),
            chain: None,
            rel: None,
            chain_id: None,
        }),
        P::Pdb(PdbSeqId {
            mol: b("1abc"),
            chain: Some(66),
            rel: Some(Date::Str(b("1999"))),
            chain_id: Some(b("A")),
        }),
        P::Pdb(PdbSeqId {
            mol: b("1abc"),
            chain: Some(32),
            rel: Some(Date::Std(DateStd {
                year: 2001,
                month: None,
                day: None,
                season: None,
                hour: None,
                minute: None,
                second: None,
            })),
            chain_id: None,
        }),
        P::Pdb(PdbSeqId {
            mol: b("1abc"),
            chain: Some(-191),
            rel: Some(Date::Std(DateStd {
                year: 2002,
                month: Some(12),
                day: Some(31),
                season: Some(b("spring")),
                hour: Some(23),
                minute: Some(59),
                second: Some(58),
            })),
            chain_id: None,
        }),
    ];
    for make in [
        P::Genbank as fn(TextseqId) -> P,
        P::Embl,
        P::Pir,
        P::Swissprot,
        P::Other,
        P::Ddbj,
        P::Prf,
        P::Tpg,
        P::Tpe,
        P::Tpd,
        P::Gpipe,
        P::NamedAnnotTrack,
    ] {
        ids.push(make(full()));
        ids.push(make(least()));
        ids.push(make(text(None, Some("x"), None, Some(0))));
        ids.push(make(text(Some("n"), None, Some("reviewed"), None)));
    }
    ids
}

#[test]
fn the_mask_reader_and_the_database_decoder_build_the_same_seq_id() {
    let mut compared = 0;
    for indefinite in [true, false] {
        let writer = Writer { indefinite };
        for id in seq_ids() {
            let data = writer.seq_id(&id);
            let from_mask = seq_id_from_binary(&data)
                .unwrap_or_else(|fail| panic!("{id:?} (indefinite {indefinite}): {fail:?}"));
            let mut replaced = 0;
            let mut ber = Reader::new(&data);
            let from_db = decode_seq_id(&mut ber, &mut replaced)
                .unwrap_or_else(|error| panic!("{id:?} (indefinite {indefinite}): {error}"));
            assert!(ber.at_record_end(), "{id:?}");
            assert_eq!(replaced, 0, "{id:?}");
            assert_eq!(from_mask, from_db, "{id:?} (indefinite {indefinite})");
            assert_eq!(from_mask, id, "{id:?} (indefinite {indefinite})");
            compared += 1;
        }
    }
    // Every alternative of the CHOICE is among the ids.
    let alternatives: std::collections::BTreeSet<u8> = seq_ids()
        .iter()
        .map(|id| Writer { indefinite: false }.seq_id(id)[0] & 0x1f)
        .collect();
    assert_eq!(alternatives.len(), 20);
    assert_eq!(compared, 2 * seq_ids().len());
}
