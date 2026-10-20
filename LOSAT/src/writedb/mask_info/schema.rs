//! The ASN.1 types of a mask file, as NCBI's serial library knows them: the members of each
//! SEQUENCE and the variants of each CHOICE in their declaration order (for the binary
//! encoding NCBI numbers them `[0]`, `[1]`, ... in that order, "automatic tags"), the
//! ENUMERATED values and which members may be left out.
//!
//! Only the types that a `Blast-db-mask-info` reaches are described: `Blast-mask-list`, `Seq-loc`
//! with everything below it (`Seq-id` and the types of its variants, `Int-fuzz`, `Date`, ...).

/// How the C++ class of a type stores a number: this decides the range that the reader accepts.
#[derive(Clone, Copy)]
pub(super) enum Ty {
    /// `INTEGER` stored as `int` (`Int4`).
    Int32,
    /// `INTEGER` stored as `Int8` (`Seq-id.gi`).
    Int64,
    /// `INTEGER` stored as `TSeqPos` (`Uint4`).
    Uint32,
    Bool,
    /// `VisibleString`.
    Str,
    Null,
    Enum(&'static EnumDef),
    Class(&'static ClassDef),
    Choice(&'static ChoiceDef),
    /// `SEQUENCE OF`.
    SeqOf(&'static Ty),
    /// `SET OF`.
    SetOf(&'static Ty),
}

pub(super) struct EnumDef {
    pub(super) values: &'static [(&'static str, i32)],
    /// `INTEGER { names }` (any integer is accepted) as opposed to `ENUMERATED`.
    pub(super) integer: bool,
}

pub(super) struct Member {
    pub(super) name: &'static str,
    pub(super) ty: &'static Ty,
    /// `OPTIONAL` or `DEFAULT`: the reader does not insist that the member is present.
    pub(super) optional: bool,
}

pub(super) struct ClassDef {
    /// The type name that a text file states in front of `::=`.
    pub(super) name: &'static str,
    pub(super) members: &'static [Member],
}

pub(super) struct Variant {
    pub(super) name: &'static str,
    pub(super) ty: &'static Ty,
}

pub(super) struct ChoiceDef {
    pub(super) variants: &'static [Variant],
}

const fn req(name: &'static str, ty: &'static Ty) -> Member {
    Member {
        name,
        ty,
        optional: false,
    }
}

const fn opt(name: &'static str, ty: &'static Ty) -> Member {
    Member {
        name,
        ty,
        optional: true,
    }
}

const fn var(name: &'static str, ty: &'static Ty) -> Variant {
    Variant { name, ty }
}

// NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.asn:51-79
// ```c++
// Blast-filter-program ::= INTEGER {
//     not-set             (0),
//     dust                (10),
//     seg                 (20),
//     windowmasker        (30),
//     repeat              (40),
//     other               (100),
//     max                 (255)
// }
//
// Blast-mask-list ::= SEQUENCE {
//     -- masks for a single sequence should be grouped in a Packed-seqint
//     masks SEQUENCE OF Seq-loc,
//     -- as of 01/21/2010, this field is set to false always, indicating that the
//     -- entire object (set of Seq-loc) is contained in this object
//     more BOOLEAN
// }
//
// -- Defines the masking information for a set of sequences
// Blast-db-mask-info ::= SEQUENCE {
//     algo-id      INTEGER,
//     algo-program Blast-filter-program,
//     algo-options VisibleString,
//     -- This object was originally created to allow pagination of the sequence
//     -- masks, but this feature was discontinued in 01/21/2010
//     masks        Blast-mask-list
// }
// ```
pub(super) static BLAST_FILTER_PROGRAM: EnumDef = EnumDef {
    values: &[
        ("not-set", 0),
        ("dust", 10),
        ("seg", 20),
        ("windowmasker", 30),
        ("repeat", 40),
        ("other", 100),
        ("max", 255),
    ],
    integer: true,
};
static T_FILTER_PROGRAM: Ty = Ty::Enum(&BLAST_FILTER_PROGRAM);

pub(super) static BLAST_DB_MASK_INFO: ClassDef = ClassDef {
    name: "Blast-db-mask-info",
    members: &[
        req("algo-id", &T_INT32),
        req("algo-program", &T_FILTER_PROGRAM),
        req("algo-options", &T_STR),
        req("masks", &T_BLAST_MASK_LIST),
    ],
};
static T_BLAST_MASK_LIST: Ty = Ty::Class(&BLAST_MASK_LIST);

pub(super) static BLAST_MASK_LIST: ClassDef = ClassDef {
    name: "Blast-mask-list",
    members: &[req("masks", &T_SEQ_LOC_LIST), req("more", &T_BOOL)],
};
static T_SEQ_LOC_LIST: Ty = Ty::SeqOf(&T_SEQ_LOC);

pub(super) static T_INT32: Ty = Ty::Int32;
pub(super) static T_INT64: Ty = Ty::Int64;
pub(super) static T_UINT32: Ty = Ty::Uint32;
pub(super) static T_BOOL: Ty = Ty::Bool;
pub(super) static T_STR: Ty = Ty::Str;
static T_NULL: Ty = Ty::Null;

// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:79-91
// ```c++
// Seq-loc ::= CHOICE {
//     null NULL ,           -- not placed
//     empty Seq-id ,        -- to NULL one Seq-id in a collection
//     whole Seq-id ,        -- whole sequence
//     int Seq-interval ,    -- from to
//     packed-int Packed-seqint ,
//     pnt Seq-point ,
//     packed-pnt Packed-seqpnt ,
//     mix Seq-loc-mix ,
//     equiv Seq-loc-equiv ,  -- equivalent sets of locations
//     bond Seq-bond ,
//     feat Feat-id }         -- indirect, through a Seq-feat
// ```
pub(super) static T_SEQ_LOC: Ty = Ty::Choice(&SEQ_LOC);
pub(super) static SEQ_LOC: ChoiceDef = ChoiceDef {
    variants: &[
        var("null", &T_NULL),
        var("empty", &T_SEQ_ID),
        var("whole", &T_SEQ_ID),
        var("int", &T_SEQ_INTERVAL),
        var("packed-int", &T_PACKED_SEQINT),
        var("pnt", &T_SEQ_POINT),
        var("packed-pnt", &T_PACKED_SEQPNT),
        var("mix", &T_SEQ_LOC_MIX),
        var("equiv", &T_SEQ_LOC_EQUIV),
        var("bond", &T_SEQ_BOND),
        var("feat", &T_FEAT_ID),
    ],
};
/// Variant numbers of `Seq-loc` that LOSAT reads into a [`super::MaskLoc`].
pub(super) const SEQ_LOC_INT: usize = 3;
pub(super) const SEQ_LOC_PACKED_INT: usize = 4;

// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:93-129
// ```c++
// Seq-interval ::= SEQUENCE {
//     from INTEGER ,
//     to INTEGER ,
//     strand Na-strand OPTIONAL ,
//     id Seq-id ,    -- WARNING: this used to be optional
//     fuzz-from Int-fuzz OPTIONAL ,
//     fuzz-to Int-fuzz OPTIONAL }
//
// Packed-seqint ::= SEQUENCE OF Seq-interval
//
// Seq-point ::= SEQUENCE {
//     point INTEGER ,
//     strand Na-strand OPTIONAL ,
//     id Seq-id ,     -- WARNING: this used to be optional
//     fuzz Int-fuzz OPTIONAL }
//
// Packed-seqpnt ::= SEQUENCE {
//     strand Na-strand OPTIONAL ,
//     id Seq-id ,
//     fuzz Int-fuzz OPTIONAL ,
//     points SEQUENCE OF INTEGER }
//
// Na-strand ::= ENUMERATED {          -- strand of nucleic acid
//     unknown (0) ,
//     plus (1) ,
//     minus (2) ,
//     both (3) ,                -- in forward orientation
//     both-rev (4) ,            -- in reverse orientation
//     other (255) }
//
// Seq-bond ::= SEQUENCE {         -- bond between residues
//     a Seq-point ,           -- connection to a least one residue
//     b Seq-point OPTIONAL }  -- other end may not be available
//
// Seq-loc-mix ::= SEQUENCE OF Seq-loc   -- this will hold anything
//
// Seq-loc-equiv ::= SET OF Seq-loc      -- for a set of equivalent locations
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.def:1-17
// ```c++
// [Seq-interval]
// from._type = TSeqPos
// to._type   = TSeqPos
//
// [Seq-point]
// point._type = TSeqPos
//
// [Packed-seqpnt]
// points._type   = vector
// points.E._type = TSeqPos
// ```
pub(super) static SEQ_INTERVAL: ClassDef = ClassDef {
    name: "Seq-interval",
    members: &[
        req("from", &T_UINT32),
        req("to", &T_UINT32),
        opt("strand", &T_NA_STRAND),
        req("id", &T_SEQ_ID),
        opt("fuzz-from", &T_INT_FUZZ),
        opt("fuzz-to", &T_INT_FUZZ),
    ],
};
pub(super) const SEQ_INTERVAL_FROM: usize = 0;
pub(super) const SEQ_INTERVAL_TO: usize = 1;
pub(super) const SEQ_INTERVAL_STRAND: usize = 2;
pub(super) const SEQ_INTERVAL_ID: usize = 3;
static T_SEQ_INTERVAL: Ty = Ty::Class(&SEQ_INTERVAL);
pub(super) static T_PACKED_SEQINT: Ty = Ty::SeqOf(&T_SEQ_INTERVAL);

static SEQ_POINT: ClassDef = ClassDef {
    name: "Seq-point",
    members: &[
        req("point", &T_UINT32),
        opt("strand", &T_NA_STRAND),
        req("id", &T_SEQ_ID),
        opt("fuzz", &T_INT_FUZZ),
    ],
};
static T_SEQ_POINT: Ty = Ty::Class(&SEQ_POINT);

static T_UINT32_LIST: Ty = Ty::SeqOf(&T_UINT32);
static PACKED_SEQPNT: ClassDef = ClassDef {
    name: "Packed-seqpnt",
    members: &[
        opt("strand", &T_NA_STRAND),
        req("id", &T_SEQ_ID),
        opt("fuzz", &T_INT_FUZZ),
        req("points", &T_UINT32_LIST),
    ],
};
static T_PACKED_SEQPNT: Ty = Ty::Class(&PACKED_SEQPNT);

pub(super) static NA_STRAND: EnumDef = EnumDef {
    values: &[
        ("unknown", 0),
        ("plus", 1),
        ("minus", 2),
        ("both", 3),
        ("both-rev", 4),
        ("other", 255),
    ],
    integer: false,
};
static T_NA_STRAND: Ty = Ty::Enum(&NA_STRAND);

static SEQ_BOND: ClassDef = ClassDef {
    name: "Seq-bond",
    members: &[req("a", &T_SEQ_POINT), opt("b", &T_SEQ_POINT)],
};
static T_SEQ_BOND: Ty = Ty::Class(&SEQ_BOND);

static T_SEQ_LOC_MIX: Ty = Ty::SeqOf(&T_SEQ_LOC);
static T_SEQ_LOC_EQUIV: Ty = Ty::SetOf(&T_SEQ_LOC);

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/seqfeat.asn:32-36
// ```c++
// Feat-id ::= CHOICE {
//     gibb INTEGER ,            -- geninfo backbone
//     giim Giimport-id ,        -- geninfo import
//     local Object-id ,         -- for local software use
//     general Dbtag }           -- for use by various databases
// ```
static FEAT_ID: ChoiceDef = ChoiceDef {
    variants: &[
        var("gibb", &T_INT32),
        var("giim", &T_GIIMPORT_ID),
        var("local", &T_OBJECT_ID),
        var("general", &T_DBTAG),
    ],
};
static T_FEAT_ID: Ty = Ty::Choice(&FEAT_ID);

// NCBI reference (598d8ae6): c++/src/objects/general/general.asn:106-117
// ```c++
// Int-fuzz ::= CHOICE {
//     p-m INTEGER ,                    -- plus or minus fixed amount
//     range SEQUENCE {                 -- max to min
//         max INTEGER ,
//         min INTEGER } ,
//     pct INTEGER ,                    -- % plus or minus (x10) 0-1000
//     lim ENUMERATED {                 -- some limit value
//         unk (0) ,                    -- unknown
//         gt (1) ,                     -- greater than
//         lt (2) ,                     -- less than
//         tr (3) ,                     -- space to right of position
//         tl (4) ,                     -- space to left of position
//         circle (5) ,                 -- artificial break at origin of circle
//         other (255) } ,              -- something else
//     alt SET OF INTEGER }             -- set of alternatives for the integer
// ```
// NCBI reference (598d8ae6): c++/src/objects/general/general.def:1-9
// ```c++
// [Int-fuzz]
// p-m._type       = TSeqPos
// range.max._type = TSeqPos
// range.min._type = TSeqPos
// alt.E._type     = TSeqPos
// ```
static INT_FUZZ_RANGE: ClassDef = ClassDef {
    name: "",
    members: &[req("max", &T_UINT32), req("min", &T_UINT32)],
};
static T_INT_FUZZ_RANGE: Ty = Ty::Class(&INT_FUZZ_RANGE);
static INT_FUZZ_LIM: EnumDef = EnumDef {
    values: &[
        ("unk", 0),
        ("gt", 1),
        ("lt", 2),
        ("tr", 3),
        ("tl", 4),
        ("circle", 5),
        ("other", 255),
    ],
    integer: false,
};
static T_INT_FUZZ_LIM: Ty = Ty::Enum(&INT_FUZZ_LIM);
static T_UINT32_SET: Ty = Ty::SetOf(&T_UINT32);
static INT_FUZZ: ChoiceDef = ChoiceDef {
    variants: &[
        var("p-m", &T_UINT32),
        var("range", &T_INT_FUZZ_RANGE),
        var("pct", &T_INT32),
        var("lim", &T_INT_FUZZ_LIM),
        var("alt", &T_UINT32_SET),
    ],
};
static T_INT_FUZZ: Ty = Ty::Choice(&INT_FUZZ);

// NCBI reference (598d8ae6): c++/src/objects/general/general.asn:54-80
// ```c++
// Date ::= CHOICE {
//     str VisibleString ,        -- for those unparsed dates
//     std Date-std }             -- use this if you can
//
// Date-std ::= SEQUENCE {        -- NOTE: this is NOT a unix tm struct
//     year INTEGER ,             -- full year (including 1900)
//     month INTEGER OPTIONAL ,   -- month (1-12)
//     day INTEGER OPTIONAL ,     -- day of month (1-31)
//     season VisibleString OPTIONAL ,  -- for "spring", "may-june", etc
//     hour INTEGER OPTIONAL ,    -- hour of day (0-23)
//     minute INTEGER OPTIONAL ,  -- minute of hour (0-59)
//     second INTEGER OPTIONAL }  -- second of minute (0-59)
//
// Dbtag ::= SEQUENCE {
//     db VisibleString ,          -- name of database or system
//     tag Object-id }         -- appropriate tag
//
// Object-id ::= CHOICE {
//     id INTEGER ,
//     str VisibleString }
// ```
static OBJECT_ID: ChoiceDef = ChoiceDef {
    variants: &[var("id", &T_INT32), var("str", &T_STR)],
};
pub(super) static T_OBJECT_ID: Ty = Ty::Choice(&OBJECT_ID);
pub(super) const OBJECT_ID_ID: usize = 0;

static DBTAG: ClassDef = ClassDef {
    name: "Dbtag",
    members: &[req("db", &T_STR), req("tag", &T_OBJECT_ID)],
};
static T_DBTAG: Ty = Ty::Class(&DBTAG);

static DATE_STD: ClassDef = ClassDef {
    name: "Date-std",
    members: &[
        req("year", &T_INT32),
        opt("month", &T_INT32),
        opt("day", &T_INT32),
        opt("season", &T_STR),
        opt("hour", &T_INT32),
        opt("minute", &T_INT32),
        opt("second", &T_INT32),
    ],
};
static T_DATE_STD: Ty = Ty::Class(&DATE_STD);
static DATE: ChoiceDef = ChoiceDef {
    variants: &[var("str", &T_STR), var("std", &T_DATE_STD)],
};
static T_DATE: Ty = Ty::Choice(&DATE);

// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:24-77
// ```c++
// Seq-id ::= CHOICE {
//     local Object-id ,            -- local use
//     gibbsq INTEGER ,             -- Geninfo backbone seqid
//     gibbmt INTEGER ,             -- Geninfo backbone moltype
//     giim Giimport-id ,           -- Geninfo import id
//     genbank Textseq-id ,
//     embl Textseq-id ,
//     pir Textseq-id ,
//     swissprot Textseq-id ,
//     patent Patent-seq-id ,
//     other Textseq-id ,           -- for historical reasons, 'other' = 'refseq'
//     general Dbtag ,              -- for other databases
//     gi INTEGER ,                 -- GenInfo Integrated Database
//     ddbj Textseq-id ,            -- DDBJ
//     prf Textseq-id ,             -- PRF SEQDB
//     pdb PDB-seq-id ,             -- PDB sequence
//     tpg Textseq-id ,             -- Third Party Annot/Seq Genbank
//     tpe Textseq-id ,             -- Third Party Annot/Seq EMBL
//     tpd Textseq-id ,             -- Third Party Annot/Seq DDBJ
//     gpipe Textseq-id ,           -- Internal NCBI genome pipeline processing ID
//     named-annot-track Textseq-id -- Internal named annotation tracking ID
// }
//
// Patent-seq-id ::= SEQUENCE {
//     seqid INTEGER ,         -- number of sequence in patent
//     cit Id-pat }           -- patent citation
//
// Textseq-id ::= SEQUENCE {
//     name VisibleString OPTIONAL ,
//     accession VisibleString OPTIONAL ,
//     release VisibleString OPTIONAL ,
//     version INTEGER OPTIONAL }
//
// Giimport-id ::= SEQUENCE {
//     id INTEGER ,                     -- the id to use here
//     db VisibleString OPTIONAL ,      -- dbase used in
//     release VisibleString OPTIONAL } -- the release
//
// PDB-seq-id ::= SEQUENCE {            --
//       mol PDB-mol-id ,                    -- the molecule name
//       chain INTEGER DEFAULT 32 ,          -- Deprecated:  'chain' can't support multiple character
//       rel Date OPTIONAL ,                 -- release date, month and year
//       chain-id VisibleString OPTIONAL }   -- chain identifier; length-independent generalization of 'chain'
//
// PDB-mol-id ::= VisibleString  -- name of mol, 4 chars
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.def:1-6
// ```c++
// [Seq-id]
// gi._type = ncbi::TGi
// gi._storage_type = ncbi::TIntId
// ```
// NCBI reference (598d8ae6): c++/src/objects/biblio/biblio.asn:119-124
// ```c++
// Id-pat ::= SEQUENCE {                   -- just to identify a patent
//     country VisibleString ,             -- Patent Document Country
//     id CHOICE {
//         number VisibleString ,          -- Patent Document Number
//         app-number VisibleString } ,    -- Patent Doc Appl Number
//     doc-type VisibleString OPTIONAL }   -- Patent Doc Type
// ```
static ID_PAT_ID: ChoiceDef = ChoiceDef {
    variants: &[var("number", &T_STR), var("app-number", &T_STR)],
};
static T_ID_PAT_ID: Ty = Ty::Choice(&ID_PAT_ID);
static ID_PAT: ClassDef = ClassDef {
    name: "Id-pat",
    members: &[
        req("country", &T_STR),
        req("id", &T_ID_PAT_ID),
        opt("doc-type", &T_STR),
    ],
};
static T_ID_PAT: Ty = Ty::Class(&ID_PAT);
static PATENT_SEQ_ID: ClassDef = ClassDef {
    name: "Patent-seq-id",
    members: &[req("seqid", &T_INT32), req("cit", &T_ID_PAT)],
};
static T_PATENT_SEQ_ID: Ty = Ty::Class(&PATENT_SEQ_ID);
static TEXTSEQ_ID: ClassDef = ClassDef {
    name: "Textseq-id",
    members: &[
        opt("name", &T_STR),
        opt("accession", &T_STR),
        opt("release", &T_STR),
        opt("version", &T_INT32),
    ],
};
static T_TEXTSEQ_ID: Ty = Ty::Class(&TEXTSEQ_ID);
static GIIMPORT_ID: ClassDef = ClassDef {
    name: "Giimport-id",
    members: &[
        req("id", &T_INT32),
        opt("db", &T_STR),
        opt("release", &T_STR),
    ],
};
static T_GIIMPORT_ID: Ty = Ty::Class(&GIIMPORT_ID);
static PDB_SEQ_ID: ClassDef = ClassDef {
    name: "PDB-seq-id",
    members: &[
        req("mol", &T_STR),
        opt("chain", &T_INT32),
        opt("rel", &T_DATE),
        opt("chain-id", &T_STR),
    ],
};
static T_PDB_SEQ_ID: Ty = Ty::Class(&PDB_SEQ_ID);

pub(super) static SEQ_ID: ChoiceDef = ChoiceDef {
    variants: &[
        var("local", &T_OBJECT_ID),
        var("gibbsq", &T_INT32),
        var("gibbmt", &T_INT32),
        var("giim", &T_GIIMPORT_ID),
        var("genbank", &T_TEXTSEQ_ID),
        var("embl", &T_TEXTSEQ_ID),
        var("pir", &T_TEXTSEQ_ID),
        var("swissprot", &T_TEXTSEQ_ID),
        var("patent", &T_PATENT_SEQ_ID),
        var("other", &T_TEXTSEQ_ID),
        var("general", &T_DBTAG),
        var("gi", &T_INT64),
        var("ddbj", &T_TEXTSEQ_ID),
        var("prf", &T_TEXTSEQ_ID),
        var("pdb", &T_PDB_SEQ_ID),
        var("tpg", &T_TEXTSEQ_ID),
        var("tpe", &T_TEXTSEQ_ID),
        var("tpd", &T_TEXTSEQ_ID),
        var("gpipe", &T_TEXTSEQ_ID),
        var("named-annot-track", &T_TEXTSEQ_ID),
    ],
};
pub(super) static T_SEQ_ID: Ty = Ty::Choice(&SEQ_ID);
pub(super) const SEQ_ID_LOCAL: usize = 0;
