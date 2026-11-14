//! NCBI's ASN.1 types as its serial library knows them (the type information that NCBI's
//! `datatool` generates from the ASN.1 modules and their `.def` files), for the types that a
//! PSSM file reaches (`PssmWithParameters`, psiblast and tblastn `-in_pssm`), and the one generic
//! reader of ASN.1 text that runs on them (`read`). This module owns the table (`table`, generated
//! by `scripts/asn_spec_gen.py` and checked against NCBI's generated `*_.cpp` by
//! `scripts/asn_spec_check.py`) and the generic value (`Value`); the PSSM's typed model
//! (`report::serial::pssm`) is built from the value at its typed boundary
//! (`blastinput::pssm_text`).
//!
//! NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1380-1397
//! ```c++
//! void CObjectIStreamAsn::ReadClassSequential(const CClassTypeInfo* classType,
//!                                             TObjectPtr classPtr)
//! {
//!     BEGIN_OBJECT_FRAME3(eFrameClass, classType, classPtr);
//!     StartBlock();
//!
//!     ReadClassSequentialContentsBegin(classType);
//!
//!     TMemberIndex index;
//!     while ( (index = BeginClassMember(classType,*pos)) != kInvalidMember ) {
//!         ReadClassSequentialContentsMember(classPtr);
//!     }
//!
//!     ReadClassSequentialContentsEnd(classPtr);
//!
//!     EndBlock();
//!     END_OBJECT_FRAME();
//! }
//! ```
//!
//! What the table holds of a type is what NCBI's text reader reads of it: the kind, the members
//! or variants in module order with their tags, OPTIONAL and DEFAULT marks and the delayed member
//! of `Auth-list` (`biblio.def`: `names._delay = 1`), the value lists, the C++ storage of each
//! INTEGER (`.def` `_type`: `TSeqPos` is read as Uint4, `TGi` as Int8, ...) and the ASN.1 name
//! (empty for an anonymous type, as `CTypeInfo::GetName` reports an internal type,
//! serial/typeinfo.cpp:249-252), which the reader's stack names in NCBI's exception reports.

mod hooks;
pub(crate) mod read;
// Generated (`scripts/asn_spec_gen.py`), checked byte for byte by `scripts/asn_spec_check.py`.
#[rustfmt::skip]
pub(crate) mod table;

#[cfg(test)]
mod read_tests;

pub(crate) use read::{read_document, Post, ReadError};

/// `Member::flags`: an OPTIONAL member (`->SetOptional()`).
pub(crate) const OPTIONAL: u8 = 1;
/// `Member::flags`: a member with a DEFAULT (`->SetDefault(...)`; NCBI reads it as optional).
pub(crate) const DEFAULT: u8 = 2;
/// `Member::flags`: a member that NCBI reads into a delay buffer (`->SetDelayBuffer(...)`):
/// skipped with the skip functions and kept as the text read (`Value::Delayed`).
pub(crate) const DELAYED: u8 = 4;

/// One ASN.1 type of the table.
#[derive(Debug)]
pub(crate) struct TypeInfo {
    /// The ASN.1 name; empty for an anonymous type (a member's own SEQUENCE, CHOICE, ENUMERATED,
    /// SEQUENCE OF or SET OF), as NCBI's internal type information reports it.
    pub(crate) name: &'static str,
    /// The module of a named type (for the writers of line RG item D3-2, which name it in XML).
    #[allow(dead_code)]
    pub(crate) module: &'static str,
    pub(crate) kind: Kind,
}

/// What a type is.
#[derive(Debug)]
pub(crate) enum Kind {
    /// SEQUENCE: its members in module order.
    Sequence(&'static [Member]),
    /// CHOICE: its variants in module order.
    Choice(&'static [Member]),
    /// SEQUENCE OF (`set` false) or SET OF (`set` true): the element type. A named one
    /// (`Seq-descr ::= SET OF Seqdesc`) is NCBI's implicit class, read as its list.
    List {
        element: u16,
        /// Read the same way; the writers (line RG item D3-2) tell them apart.
        #[allow(dead_code)]
        set: bool,
    },
    /// A named type that is another named type under its name (`Ref-ext ::= Seq-loc`).
    Alias(u16),
    /// INTEGER read as the C++ storage of its member.
    Integer(Storage),
    /// ENUMERATED (`integer` false) or INTEGER with named values (`integer` true).
    Enumerated {
        integer: bool,
        values: &'static [(&'static str, i32)],
    },
    Boolean,
    Null,
    Real,
    String(StringKind),
    Octets,
}

/// A member of a SEQUENCE or a variant of a CHOICE.
#[derive(Debug)]
pub(crate) struct Member {
    pub(crate) name: &'static str,
    /// The tag (NCBI's automatic tags: the position from 0), which `[n]` names in ASN.1 text.
    pub(crate) tag: u8,
    /// The type's index in `table::TYPES`.
    pub(crate) ty: u16,
    /// `OPTIONAL`, `DEFAULT`, `DELAYED`.
    pub(crate) flags: u8,
    /// The DEFAULT value as the module writes it (for the writers of line RG item D3-2: a DEFAULT
    /// member that was read is written, whatever its value).
    #[allow(dead_code)]
    pub(crate) default: Option<&'static str>,
}

impl Member {
    /// Whether NCBI's reader takes the member as optional (OPTIONAL or DEFAULT).
    pub(crate) fn optional(&self) -> bool {
        self.flags & (OPTIONAL | DEFAULT) != 0
    }
}

/// The C++ storage of an INTEGER (`.def` `_type`), which decides how NCBI reads it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Storage {
    /// `Int2` (`ReadInt2`: an Int4 that must fit).
    Int2,
    /// `int`, `TSignedSeqPos` (`ReadInt4`).
    Int4,
    /// `TSeqPos`, `unsigned` (`ReadUint4`: no `-`).
    Uint4,
    /// `TGi`, `TEntrezId`, `TIntId`, `Int8` (`ReadInt8`).
    Int8,
}

/// The kind of a string.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum StringKind {
    /// VisibleString: a byte outside 0x20..0x7E is read as `#` with NCBI's error message.
    Visible,
    /// UTF8String: the bytes as they are (`eFNP_Allow`).
    Utf8,
    /// StringStore: read as a VisibleString (`ReadStringStore` is `ReadString`).
    Store,
}

/// The type `ty` of the table.
pub(crate) fn type_info(ty: u16) -> &'static TypeInfo {
    &table::TYPES[usize::from(ty)]
}

/// The type `ty` with its aliases followed.
pub(crate) fn resolved(mut ty: u16) -> &'static TypeInfo {
    loop {
        match type_info(ty).kind {
            Kind::Alias(target) => ty = target,
            _ => return type_info(ty),
        }
    }
}

/// The members or variants of the SEQUENCE or CHOICE `ty`.
pub(crate) fn members(ty: u16) -> &'static [Member] {
    match resolved(ty).kind {
        Kind::Sequence(members) | Kind::Choice(members) => members,
        _ => &[],
    }
}

/// The position of the member or variant `name` of `ty` (a name the table must hold).
pub(crate) fn member_index(ty: u16, name: &str) -> usize {
    members(ty)
        .iter()
        .position(|member| member.name == name)
        .unwrap_or_else(|| panic!("{} has no member {name}", type_info(ty).name))
}

/// A value read under the table: the shape of its type, the members of a SEQUENCE by their
/// position in the module (those set, in order).
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// INTEGER of any storage.
    Int(i64),
    /// ENUMERATED or INTEGER with named values: the number (its name is the table's).
    Enum(i32),
    Bool(bool),
    Null,
    Real(f64),
    /// A string as NCBI's reader keeps it (`#` in place of the bytes it replaced).
    Str(Vec<u8>),
    Octets(Vec<u8>),
    /// A SEQUENCE: (member position, value) of the members set, in module order.
    Sequence(Vec<(u16, Value)>),
    /// A CHOICE: the variant's position and value.
    Choice(u16, Box<Value>),
    /// A SEQUENCE OF or SET OF, in the order read.
    List(Vec<Value>),
    /// A delayed member: the text NCBI's reader keeps, from its first byte that is not white
    /// space to the end of the value (`StartDelayBuffer`, objistrasn.cpp:1648-1652).
    Delayed(Vec<u8>),
    /// A Seq-id, in the Seq-id model that owns every Seq-id of LOSAT (`seqid::model`).
    SeqId(crate::seqid::ParsedSeqId),
}

impl Value {
    /// The member at position `index` of a SEQUENCE value, when set.
    pub(crate) fn member(&self, index: usize) -> Option<&Value> {
        match self {
            Value::Sequence(members) => members
                .iter()
                .find(|(position, _)| usize::from(*position) == index)
                .map(|(_, value)| value),
            _ => None,
        }
    }

    /// `member` taken out of the SEQUENCE value.
    pub(crate) fn take_member(&mut self, index: usize) -> Option<Value> {
        match self {
            Value::Sequence(members) => members
                .iter()
                .position(|(position, _)| usize::from(*position) == index)
                .map(|at| members.remove(at).1),
            _ => None,
        }
    }

    /// Sets the member at position `index` of a SEQUENCE value, in module order.
    pub(crate) fn set_member(&mut self, index: usize, value: Value) {
        if let Value::Sequence(members) = self {
            let index = index as u16;
            match members.binary_search_by_key(&index, |(position, _)| *position) {
                Ok(at) => members[at].1 = value,
                Err(at) => members.insert(at, (index, value)),
            }
        }
    }
}
