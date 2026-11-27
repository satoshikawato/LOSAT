//! Typed access to the generic reader's values (`report::serial::asn_types::Value`) by the
//! members and variants of NCBI's ASN.1 modules: what NCBI's C++ getters (`GetLocation()`,
//! `GetSegs().GetDenseg()`, ...) read of a value (the query's annotations of `-in_pssm`, the
//! Blast4-archive of `blast_formatter -archive`, the objects of makeblastdb's ASN.1 input).
//!
//! No NCBI counterpart of its own: the generic value keeps the members of a SEQUENCE by their
//! position in the module and the variant of a CHOICE by its position (`asn_types::Value`), and
//! the table generator names each member by a constant (`table::member`, `asn_types::M`).

use super::{members, resolved, Kind, Value, M};
use crate::seqid::ParsedSeqId;

/// A value of the type `ty` of the table.
#[derive(Clone, Copy, Debug)]
pub(crate) struct V<'a> {
    pub(crate) ty: u16,
    pub(crate) value: &'a Value,
}

impl<'a> V<'a> {
    pub(crate) fn new(ty: u16, value: &'a Value) -> Self {
        V { ty, value }
    }

    /// Whether the value is an element of a list of `CRef<>`s stated as `NULL`: the empty
    /// reference that NCBI's list holds (`Value::Null`), which no getter reads.
    pub(crate) fn is_null(self) -> bool {
        matches!(self.value, Value::Null)
    }

    /// Whether the value holds, at any depth, an element of a list of `CRef<>`s stated as `NULL`
    /// (`is_null`; `Kind::List::refs`): what a reader that dereferences no such element checks
    /// before it reads (the search strategy's import, `report::archive::read::request`). Walked
    /// from a list on the heap, as values nest to any depth (`Value`).
    pub(crate) fn holds_null_element(self) -> bool {
        let mut todo = vec![self];
        while let Some(v) = todo.pop() {
            match (&resolved(v.ty).kind, v.value) {
                (&Kind::List { element, refs, .. }, Value::List(values)) => {
                    for value in values {
                        if refs && matches!(value, Value::Null) {
                            return true;
                        }
                        todo.push(V::new(element, value));
                    }
                }
                (Kind::Sequence(_), Value::Sequence(set)) => {
                    let members = members(v.ty);
                    todo.extend(set.iter().map(|(position, value)| {
                        V::new(members[usize::from(*position)].ty, value)
                    }));
                }
                (Kind::Choice(_), Value::Choice(index, value)) => {
                    todo.push(V::new(members(v.ty)[usize::from(*index)].ty, value));
                }
                _ => {}
            }
        }
        false
    }

    /// The member `m` of the SEQUENCE, when set (`IsSet<Member>()` then `Get<Member>()`).
    #[inline]
    pub(crate) fn get(self, m: M) -> Option<V<'a>> {
        self.has(m);
        self.value
            .member(m.index())
            .map(|value| V::new(m.ty, value))
    }

    /// The mandatory member `m` of the SEQUENCE (`Get<Member>()`), which the reader has set.
    #[inline]
    pub(crate) fn field(self, m: M) -> V<'a> {
        match self.get(m) {
            Some(value) => value,
            None => not_in_table("a mandatory member", self.value),
        }
    }

    /// The members `ms` of the SEQUENCE, each when set (`get`), read in one pass over the members
    /// set: what a reader of several members of one value (a feature, an interval, a Bioseq)
    /// takes instead of one scan for each.
    #[inline]
    pub(crate) fn pick<const N: usize>(self, ms: [M; N]) -> [Option<V<'a>>; N] {
        let mut out = [None; N];
        for m in ms {
            self.has(m);
        }
        if let Value::Sequence(set) = self.value {
            for (position, value) in set {
                for (slot, m) in out.iter_mut().zip(ms) {
                    if m.at == *position {
                        *slot = Some(V::new(m.ty, value));
                    }
                }
            }
        }
        out
    }

    /// Panics unless the value's type has the member `m` (`M::is_of`): a constant of another
    /// type's module.
    #[inline]
    fn has(self, m: M) {
        if self.ty != m.of && !m.is_of(self.ty) {
            not_in_table("a member of its type", &(m, resolved(self.ty).name()));
        }
    }

    /// The variant of the CHOICE: its name (`SelectionName` of `Which()`) and value.
    pub(crate) fn choice(self) -> (&'static str, V<'a>) {
        match self.value {
            Value::Choice(index, value) => {
                let variant = &members(self.ty)[usize::from(*index)];
                (variant.name(), V::new(variant.ty, value))
            }
            other => not_in_table("a CHOICE", other),
        }
    }

    /// The elements of a SEQUENCE OF or SET OF that the reader keeps as values.
    pub(crate) fn items(self) -> impl Iterator<Item = V<'a>> + 'a {
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
    pub(crate) fn item(self, index: usize) -> Option<V<'a>> {
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
    pub(crate) fn len(self) -> usize {
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
    pub(crate) fn int(self) -> i64 {
        match self.value {
            Value::Int(n) => *n,
            Value::Enum(n) => i64::from(*n),
            other => not_in_table("an INTEGER", other),
        }
    }

    /// The element `index` of a SEQUENCE OF INTEGER or of ENUMERATED (`vector<int>`,
    /// `vector<TSeqPos>`, `vector<ENa_strand>` ...), `None` past its end.
    pub(crate) fn int_at(self, index: usize) -> Option<i64> {
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
    pub(crate) fn int_in(self, index: usize) -> i64 {
        match self.int_at(index) {
            Some(value) => value,
            None => not_in_table("an element", self.value),
        }
    }

    /// The elements of a SEQUENCE OF INTEGER or of ENUMERATED in order.
    pub(crate) fn ints(self) -> impl Iterator<Item = i64> + 'a {
        (0..self.len()).map(move |index| self.int_in(index))
    }

    /// The string's bytes.
    pub(crate) fn str(self) -> &'a [u8] {
        match self.value {
            Value::Str(bytes) => bytes,
            other => not_in_table("a string", other),
        }
    }

    /// The OCTET STRING's bytes.
    pub(crate) fn octets(self) -> &'a [u8] {
        match self.value {
            Value::Octets(bytes) => bytes,
            other => not_in_table("an OCTET STRING", other),
        }
    }

    /// The Seq-id (the reader's Seq-id model, `Value::SeqId`).
    pub(crate) fn seq_id(self) -> &'a ParsedSeqId {
        match self.value {
            Value::SeqId(id) => id,
            other => not_in_table("a Seq-id", other),
        }
    }
}

/// A mandatory member of those `V::pick` read (`Get<Member>()`, as `V::field`), which the reader
/// has set.
#[inline]
pub(crate) fn mandatory(member: Option<V<'_>>) -> V<'_> {
    match member {
        Some(value) => value,
        None => not_in_table("a mandatory member", &"a member not set"),
    }
}

/// A value that the type table does not give its member (the reader's values follow the table):
/// one panic site for every getter.
#[cold]
#[inline(never)]
fn not_in_table(what: &str, value: &dyn std::fmt::Debug) -> ! {
    unreachable!("{what} of the table, not {value:?}")
}
