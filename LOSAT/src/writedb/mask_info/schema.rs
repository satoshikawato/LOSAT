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

impl Ty {
    /// The name of a class or choice type (`GetName()` of its type information, the name that a
    /// `:Type` of NCBI's pointer syntax is compared with); empty for the other types.
    pub(super) fn type_name(&self) -> &'static str {
        match self {
            Ty::Class(def) => def.name,
            Ty::Choice(def) => def.name,
            _ => "",
        }
    }
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
    /// Whether the generated class has a set flag for the member (`SET_MEMBER`, `HaveSetFlag()`):
    /// such a member is read by `ReadWithSetFlagMember`, which adds `error while reading <name>`
    /// to a `CSerialException`. The generated classes give every member one except a member that
    /// holds a class or choice by reference (`ADD_NAMED_REF_MEMBER`, `CRef<>`): `ref_req`,
    /// `ref_opt`.
    pub(super) set_flag: bool,
    /// An `OPTIONAL` member held by `CRef<>` (`ref_opt`). The generated constructor allocates the
    /// object of a mandatory `CRef<>` member (`ref_req`), which NCBI reads in place; an optional one
    /// is empty when NCBI reads it, so `CPointerTypeInfo::ReadPointer` reads it with NCBI's pointer
    /// syntax (`NULL`, `@n`, `:Type value`; `Rd::read_pointer`).
    pub(super) pointer: bool,
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
    /// The type's name (`Object type= ` of its frame), empty for a CHOICE declared inside a
    /// SEQUENCE (`Id-pat.id`).
    pub(super) name: &'static str,
    pub(super) variants: &'static [Variant],
}

const fn req(name: &'static str, ty: &'static Ty) -> Member {
    Member {
        name,
        ty,
        optional: false,
        set_flag: true,
        pointer: false,
    }
}

const fn opt(name: &'static str, ty: &'static Ty) -> Member {
    Member {
        name,
        ty,
        optional: true,
        set_flag: true,
        pointer: false,
    }
}

// NCBI reference (598d8ae6): datatool's src/objects/seqloc/Seq_interval_.cpp:110-123 (from seqloc.asn)
// (The other generated classes of these types are alike: `ADD_NAMED_REF_MEMBER` has no set flag;
// `Blast-db-mask-info.masks`, `Seq-point.id`/`fuzz`, `Packed-seqpnt.id`/`fuzz`, `Seq-bond.a`/`b`,
// `Dbtag.tag`, `Patent-seq-id.cit`, `Id-pat.id` and `PDB-seq-id.rel` are such members; checked
// with makeblastdb 2.17.0's "error while reading" lines.)
const fn ref_req(name: &'static str, ty: &'static Ty) -> Member {
    Member {
        name,
        ty,
        optional: false,
        set_flag: false,
        pointer: false,
    }
}

// (checked with makeblastdb 2.17.0: `NULL` for `Seq-interval.fuzz-from` and `Seq-point.fuzz` leaves
// them unset, `@1` there fails in `GetRegisteredObject`; `NULL` for a mandatory `Seq-interval.id`
// is read as a Seq-id, the h_pointer_* fixtures.)
const fn ref_opt(name: &'static str, ty: &'static Ty) -> Member {
    Member {
        name,
        ty,
        optional: true,
        set_flag: false,
        pointer: true,
    }
}

const fn var(name: &'static str, ty: &'static Ty) -> Variant {
    Variant { name, ty }
}

// NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.asn:51-79
pub(super) static BLAST_FILTER_PROGRAM: EnumDef = EnumDef {
    values: &crate::seqdb::mask::FILTER_PROGRAMS,
    integer: true,
};
static T_FILTER_PROGRAM: Ty = Ty::Enum(&BLAST_FILTER_PROGRAM);

pub(super) static BLAST_DB_MASK_INFO: ClassDef = ClassDef {
    name: "Blast-db-mask-info",
    members: &[
        req("algo-id", &T_INT32),
        req("algo-program", &T_FILTER_PROGRAM),
        req("algo-options", &T_STR),
        ref_req("masks", &T_BLAST_MASK_LIST),
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
pub(super) static T_SEQ_LOC: Ty = Ty::Choice(&SEQ_LOC);
pub(super) static SEQ_LOC: ChoiceDef = ChoiceDef {
    name: "Seq-loc",
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
/// Variant numbers of `Seq-loc`.
pub(super) const SEQ_LOC_NULL: usize = 0;
pub(super) const SEQ_LOC_EMPTY: usize = 1;
pub(super) const SEQ_LOC_WHOLE: usize = 2;
pub(super) const SEQ_LOC_INT: usize = 3;
pub(super) const SEQ_LOC_PACKED_INT: usize = 4;
pub(super) const SEQ_LOC_PNT: usize = 5;
pub(super) const SEQ_LOC_PACKED_PNT: usize = 6;
pub(super) const SEQ_LOC_MIX: usize = 7;
pub(super) const SEQ_LOC_EQUIV: usize = 8;
pub(super) const SEQ_LOC_BOND: usize = 9;

// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:93-129
// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.def:1-17
pub(super) static SEQ_INTERVAL: ClassDef = ClassDef {
    name: "Seq-interval",
    members: &[
        req("from", &T_UINT32),
        req("to", &T_UINT32),
        opt("strand", &T_NA_STRAND),
        ref_req("id", &T_SEQ_ID),
        ref_opt("fuzz-from", &T_INT_FUZZ),
        ref_opt("fuzz-to", &T_INT_FUZZ),
    ],
};
pub(super) const SEQ_INTERVAL_FROM: usize = 0;
pub(super) const SEQ_INTERVAL_TO: usize = 1;
pub(super) const SEQ_INTERVAL_STRAND: usize = 2;
pub(super) const SEQ_INTERVAL_ID: usize = 3;
static T_SEQ_INTERVAL: Ty = Ty::Class(&SEQ_INTERVAL);
pub(super) static T_PACKED_SEQINT: Ty = Ty::SeqOf(&T_SEQ_INTERVAL);

pub(super) static SEQ_POINT: ClassDef = ClassDef {
    name: "Seq-point",
    members: &[
        req("point", &T_UINT32),
        opt("strand", &T_NA_STRAND),
        ref_req("id", &T_SEQ_ID),
        ref_opt("fuzz", &T_INT_FUZZ),
    ],
};
static T_SEQ_POINT: Ty = Ty::Class(&SEQ_POINT);

static T_UINT32_LIST: Ty = Ty::SeqOf(&T_UINT32);
pub(super) static PACKED_SEQPNT: ClassDef = ClassDef {
    name: "Packed-seqpnt",
    members: &[
        opt("strand", &T_NA_STRAND),
        ref_req("id", &T_SEQ_ID),
        ref_opt("fuzz", &T_INT_FUZZ),
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

pub(super) static SEQ_BOND: ClassDef = ClassDef {
    name: "Seq-bond",
    members: &[ref_req("a", &T_SEQ_POINT), ref_opt("b", &T_SEQ_POINT)],
};
static T_SEQ_BOND: Ty = Ty::Class(&SEQ_BOND);

static T_SEQ_LOC_MIX: Ty = Ty::SeqOf(&T_SEQ_LOC);
static T_SEQ_LOC_EQUIV: Ty = Ty::SetOf(&T_SEQ_LOC);

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/seqfeat.asn:32-36
static FEAT_ID: ChoiceDef = ChoiceDef {
    name: "Feat-id",
    variants: &[
        var("gibb", &T_INT32),
        var("giim", &T_GIIMPORT_ID),
        var("local", &T_OBJECT_ID),
        var("general", &T_DBTAG),
    ],
};
static T_FEAT_ID: Ty = Ty::Choice(&FEAT_ID);

// NCBI reference (598d8ae6): c++/src/objects/general/general.asn:106-120
// NCBI reference (598d8ae6): c++/src/objects/general/general.def:1-9
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
    name: "Int-fuzz",
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
static OBJECT_ID: ChoiceDef = ChoiceDef {
    name: "Object-id",
    variants: &[var("id", &T_INT32), var("str", &T_STR)],
};
pub(super) static T_OBJECT_ID: Ty = Ty::Choice(&OBJECT_ID);
pub(super) const OBJECT_ID_ID: usize = 0;

pub(super) static DBTAG: ClassDef = ClassDef {
    name: "Dbtag",
    members: &[req("db", &T_STR), ref_req("tag", &T_OBJECT_ID)],
};
static T_DBTAG: Ty = Ty::Class(&DBTAG);

pub(super) static DATE_STD: ClassDef = ClassDef {
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
pub(super) static DATE: ChoiceDef = ChoiceDef {
    name: "Date",
    variants: &[var("str", &T_STR), var("std", &T_DATE_STD)],
};
static T_DATE: Ty = Ty::Choice(&DATE);

// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:24-77
// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.def:1-6
// NCBI reference (598d8ae6): c++/src/objects/biblio/biblio.asn:119-124
pub(super) static ID_PAT_ID: ChoiceDef = ChoiceDef {
    name: "",
    variants: &[var("number", &T_STR), var("app-number", &T_STR)],
};
static T_ID_PAT_ID: Ty = Ty::Choice(&ID_PAT_ID);
pub(super) static ID_PAT: ClassDef = ClassDef {
    name: "Id-pat",
    members: &[
        req("country", &T_STR),
        ref_req("id", &T_ID_PAT_ID),
        opt("doc-type", &T_STR),
    ],
};
static T_ID_PAT: Ty = Ty::Class(&ID_PAT);
pub(super) static PATENT_SEQ_ID: ClassDef = ClassDef {
    name: "Patent-seq-id",
    members: &[req("seqid", &T_INT32), ref_req("cit", &T_ID_PAT)],
};
static T_PATENT_SEQ_ID: Ty = Ty::Class(&PATENT_SEQ_ID);
pub(super) static TEXTSEQ_ID: ClassDef = ClassDef {
    name: "Textseq-id",
    members: &[
        opt("name", &T_STR),
        opt("accession", &T_STR),
        opt("release", &T_STR),
        opt("version", &T_INT32),
    ],
};
static T_TEXTSEQ_ID: Ty = Ty::Class(&TEXTSEQ_ID);
pub(super) static GIIMPORT_ID: ClassDef = ClassDef {
    name: "Giimport-id",
    members: &[
        req("id", &T_INT32),
        opt("db", &T_STR),
        opt("release", &T_STR),
    ],
};
static T_GIIMPORT_ID: Ty = Ty::Class(&GIIMPORT_ID);
pub(super) static PDB_SEQ_ID: ClassDef = ClassDef {
    name: "PDB-seq-id",
    members: &[
        req("mol", &T_STR),
        opt("chain", &T_INT32),
        ref_opt("rel", &T_DATE),
        opt("chain-id", &T_STR),
    ],
};
static T_PDB_SEQ_ID: Ty = Ty::Class(&PDB_SEQ_ID);

pub(super) static SEQ_ID: ChoiceDef = ChoiceDef {
    name: "Seq-id",
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
