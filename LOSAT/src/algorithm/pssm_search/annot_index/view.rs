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

    /// The variant of the CHOICE: its name (`SelectionName(Which())`) and value.
    pub(super) fn choice(self) -> (&'static str, V<'a>) {
        match self.value {
            Value::Choice(index, value) => {
                let variant = &members(self.ty)[usize::from(*index)];
                (variant.name(), V::new(variant.ty, value))
            }
            other => unreachable!("a CHOICE of the table, not {other:?}"),
        }
    }

    /// The elements of a SEQUENCE OF or SET OF that the reader keeps as values.
    pub(super) fn items(self) -> impl Iterator<Item = V<'a>> + 'a {
        let element = match resolved(self.ty).kind {
            Kind::List { element, .. } => element,
            ref other => unreachable!("a list of the table, not {other:?}"),
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
            ref other => unreachable!("a list of the table, not {other:?}"),
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
            other => unreachable!("a list of the table, not {other:?}"),
        }
    }

    /// The INTEGER or ENUMERATED number.
    pub(super) fn int(self) -> i64 {
        match self.value {
            Value::Int(n) => *n,
            Value::Enum(n) => i64::from(*n),
            other => unreachable!("an INTEGER of the table, not {other:?}"),
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
                other => unreachable!("an INTEGER of the table, not {other:?}"),
            }),
            other => unreachable!("a list of INTEGER of the table, not {other:?}"),
        }
    }

    /// The string's bytes.
    pub(super) fn str(self) -> &'a [u8] {
        match self.value {
            Value::Str(bytes) => bytes,
            other => unreachable!("a string of the table, not {other:?}"),
        }
    }

    /// The OCTET STRING's bytes.
    pub(super) fn octets(self) -> &'a [u8] {
        match self.value {
            Value::Octets(bytes) => bytes,
            other => unreachable!("an OCTET STRING of the table, not {other:?}"),
        }
    }

    /// The Seq-id (the reader's Seq-id model, `Value::SeqId`).
    pub(super) fn seq_id(self) -> &'a ParsedSeqId {
        match self.value {
            Value::SeqId(id) => id,
            other => unreachable!("a Seq-id of the table, not {other:?}"),
        }
    }
}
