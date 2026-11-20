//! The request of an archive object (`Blast4-request`): what NCBI's archive reader takes from it
//! (`CRemoteBlast::x_GetRequestInfoFromFile` through `CImportStrategy`) and what blast_formatter
//! reads of its queries and subjects. The sequences, locations and PSSM stay the values read
//! (`asn_types::Value`), with typed access to what the readers ask of them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1717-1741
//! ```c++
//!     if (m_Archive->CanGetRequest())
//!     {
//!         CRef<objects::CBlast4_request> request(&m_Archive->SetRequest());
//!         CImportStrategy strategy(request);
//!         m_Program   = strategy.GetProgram();
//!         m_Service   = strategy.GetService();
//!         m_CreatedBy = strategy.GetCreatedBy();
//!         m_Queries    = strategy.GetQueries();
//!         m_AlgoOpts.Reset( strategy.GetAlgoOptions() );
//!         m_ProgramOpts.Reset( strategy.GetProgramOptions() );
//!
//!         if (strategy.GetSubject()->IsDatabase())
//!             x_SetDatabase(strategy.GetSubject()->GetDatabase());
//!         else
//!             m_SubjectSequences = strategy.GetSubject()->SetSequences();
//!
//!         if(m_Service == "psi")
//!         {
//!         	// Would have errored out in CImportStrategy if we can't get queue search
//!         	m_FormatOpts.Reset( strategy.GetWebFormatOptions());
//!         }
//!         // Ignore return value, want side effect of setting fields.
//!         GetSearchOptions();
//!         return;
//!     }
//! ```

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
    /// `ident` (`GetCreatedBy`): the version of the program that wrote the archive.
    pub(crate) ident: Option<Vec<u8>>,
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
    SeqLocList(Vec<Value>),
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
    /// a `queue-search` stops `CImportStrategy`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2683-2693
    // ```c++
    // bool CRemoteBlast::IsErrMsgArchive(void)
    // {
    // 	if((m_Archive) != NULL && m_Archive->IsSetRequest() && m_Archive->GetRequest().CanGetBody() &&
    // 	   m_Archive->GetRequest().GetBody().IsGet_request_info() &&
    // 	   m_Archive->GetRequest().GetBody().GetGet_request_info().IsSetRequest_id() &&
    // 	   (m_Archive->GetRequest().GetBody().GetGet_request_info().GetRequest_id()  == "Error")) {
    // 		return true;
    // 	}
    // 	return false;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:62-75
    // ```c++
    //     CImportStrategy::CImportStrategy(CRef<objects::CBlast4_request> request,
    //                                      bool ignore_unsupported_options)
    //     : m_Request(request),
    //       m_IgnoreUnsupportedOptions(ignore_unsupported_options)
    // {
    //
    //     if (m_Request.Empty()) {
    //          NCBI_THROW(CBlastException, eInvalidArgument, "CBlast4_request empty");
    //     }
    //     if (m_Request->CanGetBody() && !m_Request->GetBody().IsQueue_search() ) {
    //         NCBI_THROW(CBlastException, eInvalidArgument, "No body in CBlast4_request");
    //     }
    //     m_Data.reset(new CImportStrategyData);
    // }
    // ```
    pub(crate) fn from_value(mut request: Value) -> Result<Option<Request>, Stop> {
        let ty = named::BLAST4_REQUEST;
        let ident = take(&mut request, ty, "ident").map(|mut ident| take_str(&mut ident));
        let mut body = take(&mut request, ty, "body").expect("the mandatory body");
        let (variant, mut search) = take_variant(&mut body, named::BLAST4_REQUEST_BODY);
        if variant == "get-request-info" {
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
        let ty = named::BLAST4_QUEUE_SEARCH_REQUEST;
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
            (_, mut locs) => Subject::SeqLocList(take_items(&mut locs)),
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
            ident,
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
    /// `GetParamByName` the first parameter of a name in a list ([`Request::param`]).
    pub(crate) fn params(&self) -> &[Param] {
        &self.params
    }

    /// `GetParamByName(name)` of the list `set`: its first parameter named `name`.
    // NCBI reference (598d8ae6): c++/src/objects/blast/Blast4_parameters.cpp:186-198
    // ```c++
    // CRef <CBlast4_parameter>
    // CBlast4_parameters::GetParamByName(const string name) const
    // {
    // 	CRef <CBlast4_parameter> blast4ParamFound;
    //
    // 	ITERATE(Tdata,it,this->Get()) {
    // 		if((*it)->GetName() == name) {
    // 			blast4ParamFound = (*it);
    // 			break;
    // 		}
    // 	}
    // 	return blast4ParamFound;
    // }
    // ```
    pub(crate) fn param(&self, set: ParamSet, name: &[u8]) -> Option<&Param> {
        self.params
            .iter()
            .find(|param| param.set == set && param.name == name)
    }

    /// The size of the hit list (`CBlastOptionsHandle::GetHitlistSize()` of the options that
    /// `CBlastOptionsBuilder::GetSearchOptions` builds): the last `HitlistSize` of the
    /// algorithm and program options, read through [`Request::params`], else NCBI's default,
    /// `BLAST_HITLIST_SIZE` (500, blast_options.h:160, which every options handle sets). The one
    /// option this module reads: the whole option set of a request is line F2's
    /// (`read/options.rs`), which takes this over; until then, only this parameter's own error
    /// is NCBI's (any other parameter that NCBI's builder refuses is not checked here).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:298-306
    // ```c++
    //     case 'H':
    //         if (CBlast4Field::Get(eBlastOpt_HitlistSize).Match(p)) {
    //             opts.SetHitlistSize(v.GetInteger());
    //         } else if (CBlast4Field::Get(eBlastOpt_HspRangeMax).Match(p)) {
    //             m_HspRangeMax = v.GetInteger();
    //         } else {
    //             found = false;
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:531-545
    // ```c++
    //     if (! found) {
    //         if (m_IgnoreUnsupportedOptions)
    //             return;
    //
    //         string msg = "Internal: Error processing option [";
    //         msg += nm;
    //         msg += "] type [";
    //         msg += NStr::IntToString((int) v.Which());
    //         msg += "].";
    //
    //         NCBI_THROW(CRemoteBlastException,
    //                    eServiceNotAvailable,
    //                    msg);
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/blast/names.cpp:367-373
    // ```c++
    // bool CBlast4Field::Match(const CBlast4_parameter & p) const
    // {
    //     return (p.CanGetName()        &&
    //             p.CanGetValue()       &&
    //             p.GetValue().Which() == m_Type  &&
    //             p.GetName() == m_Name);
    // }
    // ```
    // (`HitlistSize` is an `e_Integer` field, names.cpp; `Which()` is the variant's position from
    // 1. The archive reader builds without `ignore_unsupported_options`.)
    pub(crate) fn hitlist_size(&self) -> Result<i32, Stop> {
        let mut size = 500;
        for param in &self.params {
            if param.set == ParamSet::Format || param.name != b"HitlistSize" {
                continue;
            }
            let value = param.value();
            match value.choice() {
                ("integer", integer) => size = integer.int() as i32,
                _ => {
                    let which = match value.value {
                        Value::Choice(index, _) => usize::from(*index) + 1,
                        _ => 0,
                    };
                    return Err(Stop::Exception(NcbiException::new(
                        ThrowSite {
                            file: "src/algo/blast/api/blast_options_builder.cpp",
                            line: 541,
                            kind: "CRemoteBlastException::eServiceNotAvailable",
                            function: Some(
                                "BLAST::ncbi::blast::CBlastOptionsBuilder::x_ProcessOneOption",
                            ),
                        },
                        format!("Internal: Error processing option [HitlistSize] type [{which}]."),
                    )));
                }
            }
        }
        Ok(size)
    }

    /// Whether the database of the request is protein (`x_SetDatabase`'s residue type).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1537-1546
    // ```c++
    //    EBlast4_residue_type rtype(eBlast4_residue_type_unknown);
    //
    //     if (m_Program == "blastp" ||
    //         m_Program == "blastx" ||
    //         (m_Program == "tblastn" && m_Service == "rpsblast")) {
    //
    //         rtype = eBlast4_residue_type_protein;
    //     } else {
    //         rtype = eBlast4_residue_type_nucleotide;
    //     }
    // ```
    pub(crate) fn db_is_protein(&self) -> bool {
        let program = self.program.as_slice();
        program == b"blastp"
            || program == b"blastx"
            || (program == b"tblastn" && self.service == b"rpsblast")
    }

    /// Whether the program translates the query (`Blast_QueryIsTranslated` of
    /// `NetworkProgram2BlastProgramType(program, service)`), which decides the frames of the
    /// masks; NCBI aborts on the `rpsblast` service with a program other than `blastp`,
    /// `tblastn` and `blastx`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1933-1962
    // ```c++
    //     EBlastProgramType retval = eBlastTypeUndefined;
    //     Int2 rv = BlastProgram2Number(program.c_str(), &retval);
    //     ...
    //     if (service == "rpsblast") {
    //
    //         if (program == "blastp") {
    //             retval = eBlastTypeRpsBlast;
    //         } else if (program == "tblastn" || program == "blastx") {
    //             retval = eBlastTypeRpsTblastn;
    //         } else {
    //             abort();
    //         }
    //
    //     }
    //
    //     if (service == "psi") {
    //         _ASSERT(program == "blastp");
    //         retval = eBlastTypePsiBlast;
    //     }
    // ```
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
    // ```c++
    //     switch (Which()) {
    //     case e_Pssm:
    //         retval = 1;
    //         break;
    //     case e_Seq_loc_list:
    //         retval = GetSeq_loc_list().size();
    //         break;
    //     case e_Bioseq_set:
    //         retval = GetBioseq_set().GetSeq_set().size();
    //         break;
    // ```
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
    // ```c++
    //     if (queries->IsPssm()) {
    //         const CSeq_entry& seq_entry = queries->GetPssm().GetQuery();
    //         if (seq_entry.IsSeq()) {
    //             query_ids.push_back(FindBestChoice(seq_entry.GetSeq().GetId(),
    //                                                CSeq_id::BestRank));
    //         } else {
    //             _ASSERT(seq_entry.IsSet());
    //             s_ExtractQueryIdsFromBioseqSet(seq_entry.GetSet(), query_ids);
    //         }
    //     } else if (queries->IsSeq_loc_list()) {
    //         query_ids.reserve(queries->GetSeq_loc_list().size());
    //         ITERATE(CBlast4_queries::TSeq_loc_list, i, queries->GetSeq_loc_list()) {
    //             CConstRef<CSeq_id> id((*i)->GetId());
    //             query_ids.push_back(id);
    //         }
    //     } else {
    //         _ASSERT(queries->IsBioseq_set());
    //         s_ExtractQueryIdsFromBioseqSet(queries->GetBioseq_set(), query_ids);
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2078-2087
    // ```c++
    // static void s_ExtractQueryIdsFromBioseqSet(const CBioseq_set& bss,
    //                                            CSearchResultSet::TQueryIdVector&
    //                                            query_ids)
    // {
    //     // sacrifice speed for protection against infinite loops
    //     CTypeConstIterator<objects::CBioseq> itr(ConstBegin(bss, eDetectLoops));
    //     for (; itr; ++itr) {
    //         query_ids.push_back(FindBestChoice(itr->GetId(), CSeq_id::BestRank));
    //     }
    // }
    // ```
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
fn entry_bioseqs<'a>(entry: V<'a>, out: &mut Vec<Bioseq<'a>>) {
    match entry.choice() {
        ("seq", bioseq) => out.push(Bioseq(bioseq)),
        (_, set) => set_bioseqs(set, out),
    }
}

/// A Bioseq of a request, as read.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Bioseq<'a>(V<'a>);

impl<'a> Bioseq<'a> {
    /// The Bioseq value.
    pub(crate) fn value(self) -> V<'a> {
        self.0
    }

    /// Its ids (`GetId()`), in order.
    pub(crate) fn ids(self) -> impl Iterator<Item = &'a ParsedSeqId> {
        self.0.field("id").items().map(|id| id.seq_id())
    }

    /// `FindBestChoice(GetId(), CSeq_id::BestRank)`: the first id of the best rank, `None` when
    /// every id ranks `kMax_Int` or there is none.
    pub(crate) fn best_id(self) -> Option<&'a ParsedSeqId> {
        find_best(self.ids(), best_rank_score)
    }

    /// The descriptors (`descr`) in order.
    fn descriptors(self) -> impl Iterator<Item = (&'static str, V<'a>)> {
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
    pub(crate) fn seq_data(self) -> Option<V<'a>> {
        self.0.field("inst").get("seq-data")
    }
}
