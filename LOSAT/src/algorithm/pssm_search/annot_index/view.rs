//! Typed access to the generic reader's values (`report::serial::asn_types::Value`) by the
//! member and variant names of NCBI's ASN.1 modules: what the object manager's C++ getters
//! (`GetLocation()`, `GetSegs().GetDenseg()`, ...) read of the query's annotations.
//!
//! No NCBI counterpart of its own: the generic value keeps the members of a SEQUENCE by their
//! position in the module and the variant of a CHOICE by its position (`asn_types::Value`), and
//! the type table names them (`asn_types::members`).

use crate::report::serial::asn_types::{member_index, members, resolved, Kind, Value};
use crate::seqid::ParsedSeqId;

/// A value of the type `ty` of the table.
#[derive(Clone, Copy, Debug)]
pub(super) struct V<'a> {
    pub(super) ty: u16,
    pub(super) value: &'a Value,
}

impl<'a> V<'a> {
    pub(super) fn new(ty: u16, value: &'a Value) -> Self {
        V { ty, value }
    }

    /// The member `name` of the SEQUENCE, when set (`IsSet<Member>()` then `Get<Member>()`).
    pub(super) fn get(self, name: &str) -> Option<V<'a>> {
        let index = member_index(self.ty, name);
        let ty = members(self.ty)[index].ty;
        self.value.member(index).map(|value| V::new(ty, value))
    }

    /// The mandatory member `name` of the SEQUENCE (`Get<Member>()`), which the reader has set.
    pub(super) fn field(self, name: &str) -> V<'a> {
        match self.get(name) {
            Some(value) => value,
            None => not_in_table("a mandatory member", self.value),
        }
    }

    /// The variant of the CHOICE: its name (`SelectionName` of `Which()`) and value.
    pub(super) fn choice(self) -> (&'static str, V<'a>) {
        match self.value {
            Value::Choice(index, value) => {
                let variant = &members(self.ty)[usize::from(*index)];
                (variant.name(), V::new(variant.ty, value))
            }
            other => not_in_table("a CHOICE", other),
        }
    }

    /// The elements of a SEQUENCE OF or SET OF that the reader keeps as values.
    pub(super) fn items(self) -> impl Iterator<Item = V<'a>> + 'a {
        let element = match resolved(self.ty).kind {
            Kind::List { element, .. } => element,
            ref other => not_in_table("a list", other),
        };
        let values: &'a [Value] = match self.value {
            Value::List(values) => values,
            _ => &[],
        };
        values.iter().map(move |value| V::new(element, value))
    }

    /// The element `index` of a SEQUENCE OF or SET OF kept as values (`operator[]`), `None` past
    /// its end.
    pub(super) fn item(self, index: usize) -> Option<V<'a>> {
        let element = match resolved(self.ty).kind {
            Kind::List { element, .. } => element,
            ref other => not_in_table("a list", other),
        };
        match self.value {
            Value::List(values) => values.get(index).map(|value| V::new(element, value)),
            _ => None,
        }
    }

    /// The number of elements of a list of any kind (`size()`).
    pub(super) fn len(self) -> usize {
        match self.value {
            Value::List(values) => values.len(),
            Value::Int4s(values) => values.len(),
            Value::Ints(values) => values.len(),
            Value::Reals(values) => values.len(),
            Value::Octets(bytes) => bytes.len(),
            other => not_in_table("a list", other),
        }
    }

    /// The INTEGER or ENUMERATED number.
    pub(super) fn int(self) -> i64 {
        match self.value {
            Value::Int(n) => *n,
            Value::Enum(n) => i64::from(*n),
            other => not_in_table("an INTEGER", other),
        }
    }

    /// The element `index` of a SEQUENCE OF INTEGER or of ENUMERATED (`vector<int>`,
    /// `vector<TSeqPos>`, `vector<ENa_strand>` ...), `None` past its end.
    pub(super) fn int_at(self, index: usize) -> Option<i64> {
        match self.value {
            Value::Int4s(values) => values.get(index).map(|&v| i64::from(v)),
            Value::Ints(values) => values.get(index).copied(),
            Value::List(values) => values.get(index).map(|value| match value {
                Value::Int(n) => *n,
                Value::Enum(n) => i64::from(*n),
                other => not_in_table("an INTEGER", other),
            }),
            other => not_in_table("a list of INTEGER", other),
        }
    }

    /// The element `index` of a SEQUENCE OF INTEGER or of ENUMERATED that the caller's sizes
    /// hold.
    pub(super) fn int_in(self, index: usize) -> i64 {
        match self.int_at(index) {
            Some(value) => value,
            None => not_in_table("an element", self.value),
        }
    }

    /// The elements of a SEQUENCE OF INTEGER or of ENUMERATED in order.
    pub(super) fn ints(self) -> impl Iterator<Item = i64> + 'a {
        (0..self.len()).map(move |index| self.int_in(index))
    }

    /// The string's bytes.
    pub(super) fn str(self) -> &'a [u8] {
        match self.value {
            Value::Str(bytes) => bytes,
            other => not_in_table("a string", other),
        }
    }

    /// The OCTET STRING's bytes.
    pub(super) fn octets(self) -> &'a [u8] {
        match self.value {
            Value::Octets(bytes) => bytes,
            other => not_in_table("an OCTET STRING", other),
        }
    }

    /// The Seq-id (the reader's Seq-id model, `Value::SeqId`).
    pub(super) fn seq_id(self) -> &'a ParsedSeqId {
        match self.value {
            Value::SeqId(id) => id,
            other => not_in_table("a Seq-id", other),
        }
    }
}

/// A value that the type table does not give its member (the reader's values follow the table):
/// one panic site for every getter.
#[cold]
#[inline(never)]
fn not_in_table(what: &str, value: &dyn std::fmt::Debug) -> ! {
    unreachable!("{what} of the table, not {value:?}")
}
