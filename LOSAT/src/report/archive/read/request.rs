//! The request of a search strategy (`Blast4-request`) as `CImportStrategy` takes it. The
//! sequences, locations and PSSM stay the values read (`asn_types::Value`), with typed access to
//! what the import asks of them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:62-75

use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::{member_index, members, Value};
use crate::seqid::ParsedSeqId;

/// The request of a search.
#[derive(Debug)]
pub(crate) struct Request {
    /// `program` and `service` (`GetProgram`, `GetService`), as stored.
    pub(crate) program: Vec<u8>,
    pub(crate) service: Vec<u8>,
    pub(crate) queries: Queries,
    pub(crate) subject: Subject,
    /// The parameters of `algorithm-options`, `program-options` and `format-options` in the
    /// order of the file (`params`).
    params: Vec<Param>,
}

/// The queries of a request (`Blast4-queries`), as read.
#[derive(Debug)]
pub(crate) enum Queries {
    /// `pssm`: the `PssmWithParameters` (not imported).
    Pssm,
    /// `seq-loc-list`: the `Seq-loc`s (not imported).
    SeqLocList,
    /// `bioseq-set`: the `Bioseq-set`.
    BioseqSet(Value),
}

/// The subject of a request (`Blast4-subject`), as read.
#[derive(Debug)]
pub(crate) enum Subject {
    /// `database`: the name as stored (`x_SetDatabase`).
    Database(Vec<u8>),
    /// `sequences`: the `Bioseq`s (`m_SubjectSequences`).
    Sequences(Vec<Value>),
    /// `seq-loc-list`: the `Seq-loc`s (not imported).
    SeqLocList,
}

/// The list of options a parameter comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ParamSet {
    /// `algorithm-options` (`GetAlgoOptions`).
    Algorithm,
    /// `program-options` (`GetProgramOptions`).
    Program,
    /// `format-options` (`GetWebFormatOptions`).
    Format,
}

/// A `Blast4-parameter`, as read.
#[derive(Debug)]
pub(crate) struct Param {
    pub(crate) set: ParamSet,
    pub(crate) name: Vec<u8>,
    /// The `Blast4-value`.
    value: Value,
}

impl Param {
    /// The `Blast4-value` (a CHOICE: `choice()` names its variant).
    pub(crate) fn value(&self) -> V<'_> {
        V::new(named::BLAST4_VALUE, &self.value)
    }
}

/// The member `name` of the SEQUENCE `value` of the type `ty`, moved out.
fn take(value: &mut Value, ty: u16, name: &str) -> Option<Value> {
    value.take_member(member_index(ty, name))
}

/// The variant of the CHOICE `value` of the type `ty`: its name and value, moved out.
fn take_variant(value: &mut Value, ty: u16) -> (&'static str, Value) {
    match value {
        Value::Choice(index, variant) => (members(ty)[usize::from(*index)].name(), variant.take()),
        other => unreachable!("a CHOICE of the table, not {other:?}"),
    }
}

/// The elements of a SEQUENCE OF value, moved out.
fn take_items(value: &mut Value) -> Vec<Value> {
    match value {
        Value::List(items) => std::mem::take(items),
        other => unreachable!("a list of values of the table, not {other:?}"),
    }
}

/// The bytes of a string value, moved out.
fn take_str(value: &mut Value) -> Vec<u8> {
    match value {
        Value::Str(bytes) => std::mem::take(bytes),
        other => unreachable!("a string of the table, not {other:?}"),
    }
}

impl Request {
    /// The request of a search strategy (`-import_search_strategy`): `CImportStrategy` over the
    /// `Blast4-request` the program read. A body that is not a `queue-search` stops it (an error
    /// archive's `get-request-info` too); its import never asks for `ident` (NCBI's exports
    /// write none; `E.*.ident` of line F2's probes: ignored), which is empty when the file has
    /// none.
    pub(crate) fn from_strategy(mut request: Value) -> Result<Request, NcbiException> {
        let ty = named::BLAST4_REQUEST;
        let mut body = take(&mut request, ty, "body").expect("the mandatory body");
        let (variant, mut search) = take_variant(&mut body, named::BLAST4_REQUEST_BODY);
        if variant != "queue-search" {
            return Err(NcbiException::new(
                ThrowSite {
                    file: "src/algo/blast/api/search_strategy.cpp",
                    line: 72,
                    kind: "CBlastException::eInvalidArgument",
                    function: Some("BLAST::ncbi::blast::CImportStrategy::CImportStrategy"),
                },
                "No body in CBlast4_request",
            ));
        }
        let ty = named::BLAST4_QUEUE_SEARCH_REQUEST;
        let mut text = |name: &str| {
            let mut value = take(&mut search, ty, name).expect("a mandatory string");
            take_str(&mut value)
        };
        let (program, service) = (text("program"), text("service"));
        let mut queries = take(&mut search, ty, "queries").expect("the mandatory queries");
        let queries = match take_variant(&mut queries, named::BLAST4_QUERIES) {
            ("pssm", _) => Queries::Pssm,
            ("seq-loc-list", _) => Queries::SeqLocList,
            (_, set) => Queries::BioseqSet(set),
        };
        let mut subject = take(&mut search, ty, "subject").expect("the mandatory subject");
        let subject = match take_variant(&mut subject, named::BLAST4_SUBJECT) {
            ("database", mut name) => Subject::Database(take_str(&mut name)),
            ("sequences", mut bioseqs) => Subject::Sequences(take_items(&mut bioseqs)),
            (_, _) => Subject::SeqLocList,
        };
        let mut params = Vec::new();
        for (set, name) in [
            (ParamSet::Algorithm, "algorithm-options"),
            (ParamSet::Program, "program-options"),
            (ParamSet::Format, "format-options"),
        ] {
            let Some(mut list) = take(&mut search, ty, name) else {
                continue;
            };
            params.extend(take_items(&mut list).into_iter().map(|mut param| {
                let ty = named::BLAST4_PARAMETER;
                let name = take_str(&mut take(&mut param, ty, "name").expect("the name"));
                let value = take(&mut param, ty, "value").expect("the value");
                Param { set, name, value }
            }));
        }
        Ok(Request {
            program,
            service,
            queries,
            subject,
            params,
        })
    }

    /// The parameters of the request: `algorithm-options`, `program-options` and
    /// `format-options` in the order of the file (each list in its order, a name as often as it
    /// is stored). The one accessor of the request's options: NCBI's option builder reads the
    /// first two lists in this order, each parameter in turn (`x_ProcessOptions`).
    pub(crate) fn params(&self) -> &[Param] {
        &self.params
    }
}

/// The Bioseqs of a Bioseq-set's entries, depth first (`CTypeConstIterator<CBioseq>`: in the
/// types of the table, a Bioseq is only the `seq` of a `Seq-entry`, and the entries of a set
/// are its `seq-set`).
fn set_bioseqs<'a>(set: V<'a>, out: &mut Vec<Bioseq<'a>>) {
    let mut todo: Vec<V<'a>> = set.field("seq-set").items().collect();
    todo.reverse();
    while let Some(entry) = todo.pop() {
        match entry.choice() {
            ("seq", bioseq) => out.push(Bioseq(bioseq)),
            (_, set) => {
                let at = todo.len();
                todo.extend(set.field("seq-set").items());
                todo[at..].reverse();
            }
        }
    }
}

/// The Bioseqs of a Seq-entry (`entry_bioseqs` of its set, or its Bioseq).
pub(super) fn entry_bioseqs<'a>(entry: V<'a>, out: &mut Vec<Bioseq<'a>>) {
    match entry.choice() {
        ("seq", bioseq) => out.push(Bioseq(bioseq)),
        (_, set) => set_bioseqs(set, out),
    }
}

/// A Bioseq of a request, as read.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Bioseq<'a>(V<'a>);

impl<'a> Bioseq<'a> {
    /// The Bioseq of a value of type `Bioseq` (a subject of a request).
    pub(super) fn of(value: V<'a>) -> Self {
        Bioseq(value)
    }

    /// The Bioseq value.
    pub(crate) fn value(self) -> V<'a> {
        self.0
    }

    /// Its ids (`GetId()`), in order.
    pub(crate) fn ids(self) -> impl Iterator<Item = &'a ParsedSeqId> {
        self.0.field("id").items().map(|id| id.seq_id())
    }

    /// The descriptors (`descr`) in order.
    pub(super) fn descriptors(self) -> impl Iterator<Item = (&'static str, V<'a>)> {
        self.0
            .get("descr")
            .into_iter()
            .flat_map(|descr| descr.items())
            .map(|desc| desc.choice())
    }

    /// The first `title` descriptor, as stored.
    pub(crate) fn title(self) -> Option<&'a [u8]> {
        self.descriptors()
            .find(|(variant, _)| *variant == "title")
            .map(|(_, title)| title.str())
    }

    /// `inst.mol` (`Seq-inst.mol`, ENUMERATED: 0 not-set, 1 dna, 2 rna, 3 aa, 4 na, 255 other).
    pub(crate) fn molecule(self) -> i64 {
        self.0.field("inst").field("mol").int()
    }

    /// `inst.length`, when set.
    pub(crate) fn length(self) -> Option<i64> {
        self.0
            .field("inst")
            .get("length")
            .map(|length| length.int())
    }
}
