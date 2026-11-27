//! The request of an archive object (`Blast4-request`): what NCBI's archive reader takes from it
//! (`CRemoteBlast::x_GetRequestInfoFromFile` through `CImportStrategy`) and what blast_formatter
//! reads of its queries and subjects. The sequences, locations and PSSM stay the values read
//! (`asn_types::Value`), with typed access to what the readers ask of them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1717-1741

use super::Stop;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::{member_index, members, Value};
use crate::seqid::rank::{best_rank_score, find_best};
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
    /// `pssm`: the `PssmWithParameters`.
    Pssm(Value),
    /// `seq-loc-list`: the `Seq-loc`s.
    SeqLocList(Vec<Value>),
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
    /// `seq-loc-list`: the `Seq-loc`s. NCBI's archive reader does not read them: its
    /// `SetSequences()` makes the subject an empty list of sequences, in the request it keeps
    /// too, so blast_formatter has no subject.
    SeqLocList,
}

/// The NCBI code that reads a request: the archive reader (`CRemoteBlast::LoadFromArchive`,
/// then `x_GetRequestInfoFromFile`) or the strategy import (`CImportStrategy`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Reader {
    Archive,
    Strategy,
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
    /// The request of a `Blast4-request` value: `None` for an error archive's
    /// (`IsErrMsgArchive`: a `get-request-info` whose `request-id` is "Error"); a body that is not
    /// a `queue-search` stops `CImportStrategy`, and then a request without `ident`
    /// `GetCreatedBy` (`x_GetRequestInfoFromFile` asks for the program and service, then for it).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2683-2693
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:62-75
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:164-168
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:156-158
    // (`GetProgram` and `GetService` read the mandatory members `program` and `service` of the
    // `queue-search`.)
    pub(crate) fn from_value(request: Value) -> Result<Option<Request>, Stop> {
        Request::read(request, Reader::Archive)
    }

    /// The request of a search strategy (`-import_search_strategy`): `CImportStrategy` over the
    /// `Blast4-request` the program read. A body that is not a `queue-search` stops it (an error
    /// archive's `get-request-info` too); its import never asks for `ident` (NCBI's exports
    /// write none; `E.*.ident` of line F2's probes: ignored). A search with an element stated as
    /// `NULL` in a list is refused (`Stop::NullElement`).
    pub(crate) fn from_strategy(request: Value) -> Result<Request, Stop> {
        match Request::read(request, Reader::Strategy) {
            Ok(request) => Ok(request.expect("a strategy's request is a search")),
            Err(stop @ (Stop::Exception(_) | Stop::NullElement)) => Err(stop),
            Err(stop) => unreachable!(
                "a strategy's import stops only with an exception or a NULL element: {stop:?}"
            ),
        }
    }

    fn read(mut request: Value, reader: Reader) -> Result<Option<Request>, Stop> {
        let ty = named::BLAST4_REQUEST;
        // `ident` (`GetCreatedBy`, the version of the program that wrote the archive) changes no
        // output (probes P8.001-012); NCBI's archive reader stops on a search request without it.
        let has_ident = take(&mut request, ty, "ident").is_some();
        let mut body = take(&mut request, ty, "body").expect("the mandatory body");
        let (variant, mut search) = take_variant(&mut body, named::BLAST4_REQUEST_BODY);
        if reader == Reader::Archive && variant == "get-request-info" {
            let id = V::new(named::BLAST4_GET_REQUEST_INFO_REQUEST, &search).field("request-id");
            if id.str() == b"Error" {
                return Ok(None);
            }
        }
        if variant != "queue-search" {
            return Err(Stop::Exception(NcbiException::new(
                ThrowSite {
                    file: "src/algo/blast/api/search_strategy.cpp",
                    line: 72,
                    kind: "CBlastException::eInvalidArgument",
                    function: Some("BLAST::ncbi::blast::CImportStrategy::CImportStrategy"),
                },
                "No body in CBlast4_request",
            )));
        }
        if reader == Reader::Archive && !has_ident {
            return Err(Stop::Unassigned {
                type_name: "Blast4-request",
                member: "ident",
            });
        }
        let ty = named::BLAST4_QUEUE_SEARCH_REQUEST;
        if V::new(ty, &search).holds_null_element() {
            return Err(Stop::NullElement);
        }
        let mut text = |name: &str| {
            let mut value = take(&mut search, ty, name).expect("a mandatory string");
            take_str(&mut value)
        };
        let (program, service) = (text("program"), text("service"));
        let mut queries = take(&mut search, ty, "queries").expect("the mandatory queries");
        let queries = match take_variant(&mut queries, named::BLAST4_QUERIES) {
            ("pssm", pssm) => Queries::Pssm(pssm),
            ("seq-loc-list", mut locs) => Queries::SeqLocList(take_items(&mut locs)),
            (_, set) => Queries::BioseqSet(set),
        };
        let mut subject = take(&mut search, ty, "subject").expect("the mandatory subject");
        let subject = match take_variant(&mut subject, named::BLAST4_SUBJECT) {
            ("database", mut name) => Subject::Database(take_str(&mut name)),
            ("sequences", mut bioseqs) => Subject::Sequences(take_items(&mut bioseqs)),
            _ => Subject::SeqLocList,
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
        Ok(Some(Request {
            program,
            service,
            queries,
            subject,
            params,
        }))
    }

    /// The parameters of the request: `algorithm-options`, `program-options` and
    /// `format-options` in the order of the file (each list in its order, a name as often as it
    /// is stored). The one accessor of the request's options: NCBI's option builder reads the
    /// first two lists in this order, each parameter in turn (`x_ProcessOptions`), and
    /// `GetParamByName` the first parameter of a name in a list.
    pub(crate) fn params(&self) -> &[Param] {
        &self.params
    }

    /// `GetParamByName` of the list `set`: its first parameter named `name`.
    // NCBI reference (598d8ae6): c++/src/objects/blast/Blast4_parameters.cpp:186-198
    pub(crate) fn param(&self, set: ParamSet, name: &[u8]) -> Option<&Param> {
        self.params
            .iter()
            .find(|param| param.set == set && param.name == name)
    }

    /// Whether the database of the request is protein (`x_SetDatabase`'s residue type).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1537-1546
    pub(crate) fn db_is_protein(&self) -> bool {
        let program = self.program.as_slice();
        program == b"blastp"
            || program == b"blastx"
            || (program == b"tblastn" && self.service == b"rpsblast")
    }

    /// Whether the program translates the query (`Blast_QueryIsTranslated` of the program type that
    /// `NetworkProgram2BlastProgramType` gives the program and service), which decides the frames
    /// of the masks; NCBI aborts on the `rpsblast` service with a program other than `blastp`,
    /// `tblastn` and `blastx`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1933-1962
    // (`BlastProgram2Number`, blast_util.c:278-310, compares the names without case; the
    // programs whose type has `TRANSLATED_QUERY_MASK` are blastx, tblastx and rpstblastn,
    // blast_program.h:72-89. An unknown name is `eBlastTypeUndefined`, 0.)
    pub(crate) fn query_is_translated(&self) -> Result<bool, Stop> {
        let program = self.program.as_slice();
        if self.service == b"rpsblast" {
            return match program {
                b"blastp" => Ok(false),
                b"tblastn" | b"blastx" => Ok(true),
                _ => Err(Stop::Abort(
                    "NetworkProgram2BlastProgramType: the rpsblast service with another program",
                )),
            };
        }
        if self.service == b"psi" {
            return Ok(false);
        }
        Ok(["blastx", "tblastx", "rpstblastn"]
            .iter()
            .any(|name| program.eq_ignore_ascii_case(name.as_bytes())))
    }
}

impl Queries {
    /// `CBlast4_queries::GetNumQueries()`: 1 for a PSSM, the locations, the entries at the top
    /// of the Bioseq-set.
    // NCBI reference (598d8ae6): c++/src/objects/blast/Blast4_queries.cpp:57-77
    pub(crate) fn num_queries(&self) -> usize {
        match self {
            Queries::Pssm(_) => 1,
            Queries::SeqLocList(locs) => locs.len(),
            Queries::BioseqSet(set) => V::new(named::BIOSEQ_SET, set).field("seq-set").len(),
        }
    }

    /// The PSSM of a `pssm` query.
    pub(crate) fn pssm(&self) -> Option<V<'_>> {
        match self {
            Queries::Pssm(pssm) => Some(V::new(named::PSSMWITHPARAMETERS, pssm)),
            _ => None,
        }
    }

    /// The locations of a `seq-loc-list` query.
    pub(crate) fn seq_locs(&self) -> impl Iterator<Item = V<'_>> {
        let locs: &[Value] = match self {
            Queries::SeqLocList(locs) => locs,
            _ => &[],
        };
        locs.iter().map(|loc| V::new(named::SEQ_LOC, loc))
    }

    /// The Bioseqs of the queries in the order of NCBI's iterator (`CTypeConstIterator<CBioseq>`
    /// over the Bioseq-set, or over the PSSM's query when it is a set): every Bioseq of the set's
    /// entries, depth first. None for a `seq-loc-list`; a PSSM without a query has none.
    pub(crate) fn bioseqs(&self) -> Vec<Bioseq<'_>> {
        let mut out = Vec::new();
        match self {
            Queries::BioseqSet(set) => set_bioseqs(V::new(named::BIOSEQ_SET, set), &mut out),
            Queries::Pssm(pssm) => {
                if let Some(entry) = V::new(named::PSSMWITHPARAMETERS, pssm)
                    .field("pssm")
                    .get("query")
                {
                    entry_bioseqs(entry, &mut out);
                }
            }
            Queries::SeqLocList(_) => {}
        }
        out
    }

    /// `x_ExtractQueryIds`: the id of each query, in order (`None`: NCBI's null `CRef`, of a
    /// Bioseq whose ids are all ranked worst or a location without one id).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2089-2117
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2078-2087
    // (A PSSM without a query: `GetQuery()` of an unset member.)
    pub(crate) fn query_ids(&self) -> Result<Vec<Option<&ParsedSeqId>>, Stop> {
        if let Queries::Pssm(pssm) = self {
            let query = V::new(named::PSSMWITHPARAMETERS, pssm)
                .field("pssm")
                .get("query");
            if query.is_none() {
                return Err(Stop::Unassigned {
                    type_name: "Pssm",
                    member: "query",
                });
            }
        }
        if let Queries::SeqLocList(_) = self {
            return Ok(self.seq_locs().map(super::seq_loc_id).collect());
        }
        Ok(self
            .bioseqs()
            .into_iter()
            .map(|bioseq| bioseq.best_id())
            .collect())
    }
}

impl Subject {
    /// The Bioseqs of a `sequences` subject.
    pub(crate) fn bioseqs(&self) -> Vec<Bioseq<'_>> {
        match self {
            Subject::Sequences(bioseqs) => bioseqs
                .iter()
                .map(|bioseq| Bioseq(V::new(named::BIOSEQ, bioseq)))
                .collect(),
            _ => Vec::new(),
        }
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

    /// `FindBestChoice` of the ids by `CSeq_id::BestRank`: the first id of the best rank, `None`
    /// when every id ranks `kMax_Int` or there is none.
    pub(crate) fn best_id(self) -> Option<&'a ParsedSeqId> {
        find_best(self.ids(), best_rank_score)
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

    /// The `DefLine` field of the first `CFastaReader` user object, as stored (the raw defline
    /// line of the FASTA reader, `NStr::PrintableString` of it).
    #[cfg(test)] // read by the formats of later pull requests (decision C3-3)
    pub(crate) fn fasta_defline(self) -> Option<&'a [u8]> {
        self.descriptors()
            .filter(|(variant, _)| *variant == "user")
            .map(|(_, user)| user)
            .find(|user| matches!(user.field("type").choice(), ("str", kind) if kind.str() == b"CFastaReader"))
            .and_then(|user| {
                user.field("data").items().find_map(|field| {
                    match (field.field("label").choice(), field.field("data").choice()) {
                        (("str", label), ("str", text)) if label.str() == b"DefLine" => {
                            Some(text.str())
                        }
                        _ => None,
                    }
                })
            })
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

    /// `inst.seq-data` as read (a `Seq-data` CHOICE of any coding).
    #[cfg(test)] // read by the formats of later pull requests (decision C3-3)
    pub(crate) fn seq_data(self) -> Option<V<'a>> {
        self.0.field("inst").get("seq-data")
    }
}
