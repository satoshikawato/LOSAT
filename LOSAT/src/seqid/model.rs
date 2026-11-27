//! The value of an NCBI Seq-id, every alternative of the CHOICE, as the header of a BLAST
//! database stores it (`seqid::ber` decodes it). The formats that write the Seq-id object
//! (outfmt 8, 9, 11, 12) write it through `report::seqid::SeqId::Db`.
//!
//! NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:24-74
//! NCBI reference (598d8ae6): c++/src/objects/general/general.asn:54-80
//! NCBI reference (598d8ae6): c++/src/objects/biblio/biblio.asn:119-124
//!
//! An OPTIONAL member is an `Option` (`None`: not set); a member with a DEFAULT (`chain`) is
//! an `Option` too, `None` when the stored value did not hold it: NCBI keeps a set flag for
//! it, reads the default for one that was not stored, and writes it by its set flag
//! (`report::seqid` quotes the rule). The INTEGER members are NCBI's `int` (`i32`) except
//! `gi`, an `Int8` (`TGi`). Strings are the bytes after NCBI's repair of a VisibleString
//! (`seqid::ber`).

/// `Object-id`: the tag of a local id and of a `Dbtag`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ObjectId {
    /// `id INTEGER`.
    Id(i32),
    /// `str VisibleString`.
    Str(Vec<u8>),
}

/// `Dbtag`: the database name and the tag of a `general` id (`gnl|db|tag`).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Dbtag {
    pub db: Vec<u8>,
    pub tag: ObjectId,
}

/// `Textseq-id`: the accession-based ids (`gb`, `emb`, `pir`, `sp`, `ref`, `dbj`, `prf`, `tpg`,
/// `tpe`, `tpd`, `gpp`, `nat`).
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TextseqId {
    pub name: Option<Vec<u8>>,
    pub accession: Option<Vec<u8>>,
    pub release: Option<Vec<u8>>,
    pub version: Option<i32>,
}

/// `Giimport-id`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GiimportId {
    pub id: i32,
    pub db: Option<Vec<u8>>,
    pub release: Option<Vec<u8>>,
}

/// The `id` CHOICE of `Id-pat`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PatentNumber {
    /// `number VisibleString`: the patent document number.
    Number(Vec<u8>),
    /// `app-number VisibleString`: the application number.
    AppNumber(Vec<u8>),
}

/// `Id-pat`: the patent citation of a `patent` id.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct IdPat {
    pub country: Vec<u8>,
    pub id: PatentNumber,
    pub doc_type: Option<Vec<u8>>,
}

/// `Patent-seq-id`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PatentSeqId {
    /// The number of the sequence in the patent.
    pub seqid: i32,
    pub cit: IdPat,
}

/// `Date-std`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DateStd {
    pub year: i32,
    pub month: Option<i32>,
    pub day: Option<i32>,
    pub season: Option<Vec<u8>>,
    pub hour: Option<i32>,
    pub minute: Option<i32>,
    pub second: Option<i32>,
}

/// `Date`: the release date of a `pdb` id.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Date {
    Str(Vec<u8>),
    Std(DateStd),
}

/// `PDB-seq-id`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PdbSeqId {
    /// `PDB-mol-id`: the molecule name.
    pub mol: Vec<u8>,
    /// `chain INTEGER DEFAULT 32`, `None` when not stored (it reads as 32).
    pub chain: Option<i32>,
    pub rel: Option<Date>,
    pub chain_id: Option<Vec<u8>>,
}

/// The default of `PDB-seq-id.chain` (`chain INTEGER DEFAULT 32`, seqloc.asn:67, quoted
/// in the module).
pub const PDB_DEFAULT_CHAIN: i32 = 32;

impl PdbSeqId {
    /// `CPDB_seq_id::GetChain`: the stored chain, else the default (datatool's generated
    /// getter returns the member, which holds the DEFAULT until a value is read).
    pub fn chain_or_default(&self) -> i32 {
        self.chain.unwrap_or(PDB_DEFAULT_CHAIN)
    }
}

/// A Seq-id of every alternative of the CHOICE, owned: the ids that a definition line of a
/// BLAST database built with `-parse_seqids` holds (`seqid::defline_set::BlastDefline::ids`).
/// Ordered by its alternatives and members (`Ord`) for the sets that hold equal ids once
/// (psiblast's inclusion sets, `algorithm::psiblast::rounds::SetKey`).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ParsedSeqId {
    Local(ObjectId),
    Gibbsq(i32),
    Gibbmt(i32),
    Giim(GiimportId),
    Genbank(TextseqId),
    Embl(TextseqId),
    Pir(TextseqId),
    Swissprot(TextseqId),
    Patent(PatentSeqId),
    /// `other`: a RefSeq id (`ref|`).
    Other(TextseqId),
    General(Dbtag),
    Gi(i64),
    Ddbj(TextseqId),
    Prf(TextseqId),
    Pdb(PdbSeqId),
    Tpg(TextseqId),
    Tpe(TextseqId),
    Tpd(TextseqId),
    Gpipe(TextseqId),
    NamedAnnotTrack(TextseqId),
}

/// The alternatives of `Seq-id` in the order of seqloc.asn:24-45: the 0-based position is the
/// alternative's context tag in BER.
pub(crate) const SEQ_ID_ALTERNATIVES: [&str; 20] = [
    "local",
    "gibbsq",
    "gibbmt",
    "giim",
    "genbank",
    "embl",
    "pir",
    "swissprot",
    "patent",
    "other",
    "general",
    "gi",
    "ddbj",
    "prf",
    "pdb",
    "tpg",
    "tpe",
    "tpd",
    "gpipe",
    "named-annot-track",
];

impl ParsedSeqId {
    /// The 0-based position of the alternative in the CHOICE (its BER context tag).
    pub fn alternative(&self) -> u8 {
        match self {
            ParsedSeqId::Local(_) => 0,
            ParsedSeqId::Gibbsq(_) => 1,
            ParsedSeqId::Gibbmt(_) => 2,
            ParsedSeqId::Giim(_) => 3,
            ParsedSeqId::Genbank(_) => 4,
            ParsedSeqId::Embl(_) => 5,
            ParsedSeqId::Pir(_) => 6,
            ParsedSeqId::Swissprot(_) => 7,
            ParsedSeqId::Patent(_) => 8,
            ParsedSeqId::Other(_) => 9,
            ParsedSeqId::General(_) => 10,
            ParsedSeqId::Gi(_) => 11,
            ParsedSeqId::Ddbj(_) => 12,
            ParsedSeqId::Prf(_) => 13,
            ParsedSeqId::Pdb(_) => 14,
            ParsedSeqId::Tpg(_) => 15,
            ParsedSeqId::Tpe(_) => 16,
            ParsedSeqId::Tpd(_) => 17,
            ParsedSeqId::Gpipe(_) => 18,
            ParsedSeqId::NamedAnnotTrack(_) => 19,
        }
    }

    /// The ASN.1 name of the alternative (`genbank`, `other`, ...).
    pub fn alternative_name(&self) -> &'static str {
        SEQ_ID_ALTERNATIVES[usize::from(self.alternative())]
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:169-199
    /// `CSeq_id::GetTextseq_Id`: the Textseq-id of the twelve accession-based alternatives.
    pub fn textseq(&self) -> Option<&TextseqId> {
        match self {
            ParsedSeqId::Genbank(id)
            | ParsedSeqId::Embl(id)
            | ParsedSeqId::Ddbj(id)
            | ParsedSeqId::Pir(id)
            | ParsedSeqId::Swissprot(id)
            | ParsedSeqId::Other(id)
            | ParsedSeqId::Prf(id)
            | ParsedSeqId::Tpg(id)
            | ParsedSeqId::Tpe(id)
            | ParsedSeqId::Tpd(id)
            | ParsedSeqId::Gpipe(id)
            | ParsedSeqId::NamedAnnotTrack(id) => Some(id),
            _ => None,
        }
    }

    /// `textseq`, to change.
    pub fn textseq_mut(&mut self) -> Option<&mut TextseqId> {
        match self {
            ParsedSeqId::Genbank(id)
            | ParsedSeqId::Embl(id)
            | ParsedSeqId::Ddbj(id)
            | ParsedSeqId::Pir(id)
            | ParsedSeqId::Swissprot(id)
            | ParsedSeqId::Other(id)
            | ParsedSeqId::Prf(id)
            | ParsedSeqId::Tpg(id)
            | ParsedSeqId::Tpe(id)
            | ParsedSeqId::Tpd(id)
            | ParsedSeqId::Gpipe(id)
            | ParsedSeqId::NamedAnnotTrack(id) => Some(id),
            _ => None,
        }
    }

    /// The GI of a `gi` id.
    pub fn gi(&self) -> Option<i64> {
        match *self {
            ParsedSeqId::Gi(gi) => Some(gi),
            _ => None,
        }
    }
}
