//! The reader of a BLAST archive file (`-outfmt 11`, `Blast4-archive` objects one after another)
//! for `blast_formatter -archive`: the objects as NCBI's `CRemoteBlast` loads them
//! (`LoadFromArchive`), each as its stored messages, its typed request (`request`) and its
//! results (`results`).
//!
//! The objects are read by the generic reader of the type table (`report::serial::asn_types`:
//! the encoding guessed by `stream::guess_encoding`, the text read by `stream::ObjectStream`);
//! this module makes the typed model of each value. It models what NCBI's archive reader and
//! blast_formatter read of the archive; the members they ignore stay in the values (the Bioseqs,
//! Seq-locs, PSSM and Seq-aligns are kept as read, for the formats that write them again).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1657-1664
//! ```c++
//! bool
//! CRemoteBlast::LoadFromArchive()
//! {
//!       if (m_ObjectStream->EndOfData())
//!          return false;
//!
//!       m_Archive.Reset(new CBlast4_archive);
//!       *m_ObjectStream >> *m_Archive;
//! ```
//!
//! What NCBI makes of the parts of a stored request (a search strategy of
//! `-import_search_strategy`):
//!
//! - `options`: the request's options lists as the program's options (`CBlastOptionsBuilder`);
//! - `queries`: the request's query Bioseqs as the FASTA text NCBI writes and reads back
//!   (`s_ImportQueries`, `CFastaOstream`).

pub mod options;
pub mod queries;
pub(crate) mod request;
pub(crate) mod results;

#[cfg(test)]
mod read_tests;
#[cfg(test)]
mod strategy_tests;

use std::io::Read;

use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_types::stream::ObjectStream;
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::{member_index, Post, ReadError, Value};
use crate::seqid::ParsedSeqId;

use request::Request;
use results::Results;

/// Where NCBI's code stops on an archive that NCBI's own writer does not make: what the caller
/// reports in its place.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Stop {
    /// A `CException` that NCBI throws: its place, class and message.
    Exception(NcbiException),
    /// A `C<type>::Get<variant>()` of a CHOICE that holds another variant
    /// (`CInvalidChoiceSelection`, thrown from code that `datatool` generates): the module, the
    /// type, the variant asked for and the one held.
    InvalidSelection {
        module: &'static str,
        type_name: &'static str,
        required: &'static str,
        current: &'static str,
    },
    /// A `C<type>::Get<member>()` of an OPTIONAL member that is not set
    /// (`CSerialObject::ThrowUnassigned`): the type and member.
    Unassigned {
        type_name: &'static str,
        member: &'static str,
    },
    /// A `CRef` or `CConstRef` that is null is dereferenced (`CCoreException::eNullPtr`, whose
    /// report prints memory addresses): what NCBI dereferences.
    NullPointer(&'static str),
    /// NCBI's code calls `abort()`: where.
    Abort(&'static str),
    /// NCBI's code reads memory its input does not decide (an element past the end of a
    /// `vector`, `front()` of an empty `list`): where.
    Undefined(&'static str),
}

/// The severity of a stored message (`EDiagSev`, the `code` of a `Blast4-error`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Severity {
    Info,
    Warning,
    Error,
    Critical,
    /// NCBI's post of this severity ends the program.
    Fatal,
    Trace,
}

/// A stored message (`Blast4-error` with a `message`), as `LoadFromArchive` posts it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Message {
    pub(crate) severity: Severity,
    /// The text posted: the message, after "Unknown Error Code: " for a code that is no
    /// severity.
    pub(crate) text: Vec<u8>,
}

/// One object of the archive.
#[derive(Debug)]
pub(crate) struct Archive {
    /// The stored messages in order (posted as soon as the object is read, before anything else
    /// is made of it).
    pub(crate) messages: Vec<Message>,
    /// What the object holds: an error archive (`Ok(None)`), a search, or where NCBI stops after
    /// posting the messages.
    pub(crate) search: Result<Option<Search>, Stop>,
}

/// The request and results of a search.
#[derive(Debug)]
pub(crate) struct Search {
    pub(crate) request: Request,
    pub(crate) results: Results,
}

/// The objects of an archive stream (`CRemoteBlast(istr)` with its object stream), one at a
/// time.
pub(crate) struct ArchiveReader<R: Read> {
    objects: ObjectStream<R>,
}

impl<R: Read> ArchiveReader<R> {
    /// The archive objects of a text stream.
    pub(crate) fn new(objects: ObjectStream<R>) -> Self {
        ArchiveReader { objects }
    }

    /// `LoadFromArchive()`: `None` at the end of the data (`EndOfData`), else the next object
    /// read (or the reader's error, after which the stream ends) and the messages the reader
    /// posted while it read.
    pub(crate) fn next_archive(&mut self) -> Option<(Result<Archive, ReadError>, Vec<Post>)> {
        if self.objects.end_of_data() {
            return None;
        }
        let (value, posts) = self.objects.read(named::BLAST4_ARCHIVE);
        Some((value.map(Archive::from_value), posts))
    }
}

impl Archive {
    /// The archive of a `Blast4-archive` value: its messages, then what `LoadFromArchive` makes
    /// of the request (none for an error archive, `IsErrMsgArchive`) and what blast_formatter
    /// reads of the results.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1700-1705
    // ```c++
    //       if (!IsErrMsgArchive()) {
    //     	  x_GetRequestInfoFromFile(); // update info.
    //       }
    //       return true;
    // ```
    pub(crate) fn from_value(mut value: Value) -> Archive {
        let ty = named::BLAST4_ARCHIVE;
        let messages = value
            .member(member_index(ty, "messages"))
            .map(|list| messages(V::new(member_type(ty, "messages"), list)))
            .unwrap_or_default();
        let request = value
            .take_member(member_index(ty, "request"))
            .expect("the mandatory request");
        let results = value
            .take_member(member_index(ty, "results"))
            .expect("the mandatory results");
        let search = Request::from_value(request).map(|request| {
            request.map(|request| Search {
                request,
                results: Results::from_value(results),
            })
        });
        Archive { messages, search }
    }
}

/// The type of the member `name` of `ty`.
pub(crate) fn member_type(ty: u16, name: &str) -> u16 {
    crate::report::serial::asn_types::members(ty)[member_index(ty, name)].ty
}

/// The stored messages as `LoadFromArchive` posts them: each with a `message`, its `code` taken
/// as the severity.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1667-1698
// ```c++
//       if(m_Archive->IsSetMessages()) {
//     	  const list< CRef< CBlast4_error > >  & msgs = m_Archive->GetMessages();
//     	  ITERATE(list< CRef< CBlast4_error > >, itr, msgs) {
//     		  if((*itr)->IsSetMessage()) {
//     			  string msg = (*itr)->GetMessage();
//     			  if((*itr)->IsSetCode()) {
//     				  switch ((*itr)->GetCode()) {
//     			  	  case EDiagSev::eDiag_Info :
//     			  		  LOG_POST(Info << msg);
//     			  		  break;
//     			  	  ...
//     				  case EDiagSev::eDiag_Trace :
//     			  		  LOG_POST(Trace << msg);
//     					  break;
//     				  default :
//     			  		  LOG_POST(Error <<"Unknown Error Code: " << msg);
//     					  break;
//     				  }
// ```
// (`code` is mandatory, so it is always set; `EDiagSev` is Info 0, Warning 1, Error 2,
// Critical 3, Fatal 4, Trace 5, corelib/ncbidiag.hpp.)
fn messages(list: V<'_>) -> Vec<Message> {
    list.items()
        .filter_map(|error| {
            let text = error.get("message")?.str();
            let (severity, prefix): (Severity, &[u8]) = match error.field("code").int() {
                0 => (Severity::Info, b""),
                1 => (Severity::Warning, b""),
                2 => (Severity::Error, b""),
                3 => (Severity::Critical, b""),
                4 => (Severity::Fatal, b""),
                5 => (Severity::Trace, b""),
                _ => (Severity::Error, b"Unknown Error Code: "),
            };
            Some(Message {
                severity,
                text: [prefix, text].concat(),
            })
        })
        .collect()
}

/// `CSeq_loc::GetId()`: the one Seq-id of a location, `None` when it has none (`null`), when
/// its parts name ids that are not equal (`Equals`), or for a `feat`.
// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_loc.hpp:940-945
// ```c++
// inline
// const CSeq_id* CSeq_loc::GetId(void) const
// {
//     const CSeq_id* sip = NULL;
//     return CheckId(sip, false) ? sip : NULL;
// }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:3541-3630
// ```c++
// bool CSeq_loc::x_CheckId(const CSeq_id*& id, bool may_throw) const
// {
//     switch ( Which() ) {
//     case e_not_set:
//     case e_Null:
//         {
//             // no Seq-id
//             return true;
//         }
//     case e_Empty:
//         {
//             return x_UpdateId(id, &GetEmpty(), may_throw);
//         }
//     ...
//     case e_Packed_int:
//         {
//             // Check ID of each interval
//             const CPacked_seqint& ints = GetPacked_int();
//             ITERATE ( CPacked_seqint::Tdata, ii, ints.Get() ) {
//                 const CSeq_interval& loc = **ii;
//                 if ( !x_UpdateId(id, &loc.GetId(), may_throw) ) {
//                     return false;
//                 }
//             }
//             return true;
//         }
//     ...
//     case e_Mix:
//         {
//             // Check ID of each sub-location.
//             const CSeq_loc_mix& mix = GetMix();
//             ITERATE( CSeq_loc_mix::Tdata, li, mix.Get() ) {
//                 if ( !(*li)->CheckId(id, may_throw) ) {
//                     return false;
//                 }
//             }
//             return true;
//         }
//     ...
//     case e_Feat:
//     default:
//         if (may_throw) {
//             ...
//         }
//         return false;
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:403-421
// ```c++
// bool CSeq_loc::x_UpdateId(const CSeq_id*& total_id, const CSeq_id* id,
//                           bool may_throw) const
// {
//     if ( total_id == id ) {
//         return true;
//     }
//
//     if ( !total_id ) {
//         total_id = id;
//     } else if ( (id  &&  !total_id->Equals(*id)) ) {
//         if (may_throw) {
//             ...
//         } else {
//             return false;
//         }
//     }
//     return true;
// ```
// (A `mix` or `equiv` checks each part with its own id first, then that id against the total:
// the same as checking every id of the parts against one total, in any order. `bond` checks `a`,
// then `b` when set.)
pub(crate) fn seq_loc_id(loc: V<'_>) -> Option<&ParsedSeqId> {
    let mut total = None;
    let mut todo = vec![loc];
    while let Some(loc) = todo.pop() {
        let (variant, value) = loc.choice();
        let same = match variant {
            "null" => true,
            "empty" | "whole" => update_id(&mut total, value.seq_id()),
            "int" | "pnt" | "packed-pnt" => update_id(&mut total, value.field("id").seq_id()),
            "packed-int" => value
                .items()
                .all(|interval| update_id(&mut total, interval.field("id").seq_id())),
            "bond" => {
                update_id(&mut total, value.field("a").field("id").seq_id())
                    && value
                        .get("b")
                        .is_none_or(|b| update_id(&mut total, b.field("id").seq_id()))
            }
            "mix" | "equiv" => {
                // The parts in order, the first on top.
                let at = todo.len();
                todo.extend(value.items());
                todo[at..].reverse();
                true
            }
            _ => false,
        };
        if !same {
            return None;
        }
    }
    total
}

/// `x_UpdateId(total, id)` without its exception: whether `id` agrees with the id found so far
/// (the first is taken).
fn update_id<'a>(total: &mut Option<&'a ParsedSeqId>, id: &'a ParsedSeqId) -> bool {
    match total {
        None => {
            *total = Some(id);
            true
        }
        Some(total) => *total == id,
    }
}

/// `__FILE__` of `Seq_align.cpp`, `Dense_seg.cpp` and `Sparse_seg.cpp` in BLAST+'s build (the
/// `seq` library compiles them from `src/objects/seq`), as NCBI's reports name them.
const SEQ_ALIGN_CPP: &str = "src/objects/seq/../seqalign/Seq_align.cpp";
const DENSE_SEG_CPP: &str = "src/objects/seq/../seqalign/Dense_seg.cpp";
const SPARSE_SEG_CPP: &str = "src/objects/seq/../seqalign/Sparse_seg.cpp";

/// A `CSeqalignException` of `GetSeq_id` at `file`:`line`.
fn seqalign_exception(file: &'static str, line: u32, kind: &'static str, text: &str) -> Stop {
    Stop::Exception(NcbiException::new(
        ThrowSite {
            file,
            line,
            kind,
            function: None,
        },
        text.to_string(),
    ))
}

/// `CSeq_align::GetSeq_id(row)`: the Seq-id of the row `row` of an alignment, or NCBI's
/// `CSeqalignException` when the segments name none.
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:317-400
// ```c++
// const CSeq_id& CSeq_align::GetSeq_id(TDim row) const
// {
//     switch (GetSegs().Which()) {
//     case C_Segs::e_Denseg:
//         return GetSegs().GetDenseg().GetSeq_id(row);
//
//     case C_Segs::e_Sparse:
//         return GetSegs().GetSparse().GetSeq_id(row);
//
//     case C_Segs::e_Dendiag:
//         {{
//             // If segments have different number of rows, this will try
//             // to find the segment with enough rows to get the id.
//             ITERATE(CSeq_align::C_Segs::TDendiag, seg, GetSegs().GetDendiag()) {
//                 if( (*seg)->IsSetIds()  &&
//                     (size_t)row < (*seg)->GetIds().size() ) {
//                     return *(*seg)->GetIds()[row];
//                 }
//             }
//             break;
//         }}
//
//     case C_Segs::e_Std:
//         {{
//             ITERATE(CSeq_align::C_Segs::TStd, seg, GetSegs().GetStd()) {
//                 if ( (*seg)->IsSetIds()  &&
//                      (size_t)row < (*seg)->GetIds().size()) {
//                     return *((*seg)->GetIds()[row]);
//                 }
//                 else if ( (*seg)->IsSetLoc()  &&
//                      (size_t)row < (*seg)->GetLoc().size() ) {
//                     const CSeq_loc& loc = *(*seg)->GetLoc()[row];
//                     CConstRef<CSeq_id> id(loc.GetId());
//                     if (id) {
//                         return *id;
//                     }
//                 }
//             }
//         }}
//         break;
//
//     case C_Segs::e_Disc:
//         {{
//             // Try to find a sub-alignment for which we can get a
//             // Seq-id for this row.
//             ITERATE (CSeq_align_set::Tdata, sub_aln, GetSegs().GetDisc().Get()) {
//                 try {
//                     const CSeq_id& rv = (*sub_aln)->GetSeq_id(row);
//                     return rv;
//                 }
//                 catch (const CSeqalignException&) {
//                 }
//             }
//         }}
//         break;
//
//     case C_Segs::e_Spliced:
//         {{
//             const C_Segs::TSpliced& spliced_seg = GetSegs().GetSpliced();
//             if (row == 0 && spliced_seg.IsSetProduct_id()) {
//                 return spliced_seg.GetProduct_id();
//             } else if (row == 1 && spliced_seg.IsSetGenomic_id()) {
//                 return spliced_seg.GetGenomic_id();
//             }
//         }}
//         break;
//
//     default:
//         NCBI_THROW(CSeqalignException, eUnsupported,
//                    "CSeq_align::GetSeq_id() currently does not handle "
//                    "this type of alignment.");
//     }
//
//     NCBI_THROW(CSeqalignException, eInvalidRowNumber,
//                "CSeq_align::GetSeq_id(): "
//                "can not get seq-id for the row requested.");
//     // return CSeq_id();
// }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:154-162
// ```c++
// const CSeq_id& CDense_seg::GetSeq_id(TDim row) const
// {
//     if ( IsSetIds()  && (size_t)row < GetIds().size()) {
//         return *GetIds()[row];
//     }
//     NCBI_THROW(CSeqalignException, eInvalidRowNumber,
//                "CDense_seg::GetSeq_id(): "
//                "can not get seq-id for the row requested.");
// }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Sparse_seg.cpp:72-90
// ```c++
// const CSeq_id& CSparse_seg::GetSeq_id(TDim row) const
// {
//     if ( row == 0 ) {
//         // consensus sequence - first-id
//         if ( !GetRows().empty() ) {
//             const CSparse_align& align = *GetRows()[0];
//             return align.GetFirst_id();
//         }
//     }
//     else {
//         if ( size_t(row) <= GetRows().size() ) {
//             const CSparse_align& align = *GetRows()[row-1];
//             return align.GetSecond_id();
//         }
//     }
//     NCBI_THROW(CSeqalignException, eInvalidRowNumber,
//                "CSparse_seg::GetSeq_id(): "
//                "can not get seq-id for the row requested.");
// }
// ```
// (`GetSeq_id` is a `const` method, so NCBI's report names no function. The sub-alignments of
// `disc` are searched in order from a list, not by a recursion as deep as they nest: the reader
// reads any depth.)
pub(crate) fn seq_align_id(align: V<'_>, row: usize) -> Result<&ParsedSeqId, Stop> {
    // The alignments of `disc` segments still to search, the next on top; the first search is
    // the one whose exception leaves (`disc` catches those of its sub-alignments).
    let mut todo = vec![align];
    let mut outer = Ok(());
    let mut first = true;
    while let Some(align) = todo.pop() {
        match align_row_id(align, row, &mut todo) {
            Ok(Some(id)) => return Ok(id),
            Ok(None) => {}
            Err(stop) => {
                if first {
                    outer = Err(stop);
                }
            }
        }
        first = false;
    }
    match outer {
        Err(stop) => Err(stop),
        Ok(()) => Err(invalid_row()),
    }
}

/// `CSeq_align::GetSeq_id(): can not get seq-id for the row requested.` (Seq_align.cpp:396).
fn invalid_row() -> Stop {
    seqalign_exception(
        SEQ_ALIGN_CPP,
        396,
        "CSeqalignException::eInvalidRowNumber",
        "CSeq_align::GetSeq_id(): can not get seq-id for the row requested.",
    )
}

/// One alignment of `seq_align_id`: its id of the row, `None` when it has none and no exception
/// leaves it (a `disc` puts its sub-alignments on `todo`), or the exception of a segment kind
/// that throws its own.
fn align_row_id<'a>(
    align: V<'a>,
    row: usize,
    todo: &mut Vec<V<'a>>,
) -> Result<Option<&'a ParsedSeqId>, Stop> {
    let (variant, segs) = align.field("segs").choice();
    match variant {
        "denseg" => match segs.field("ids").item(row) {
            Some(id) => Ok(Some(id.seq_id())),
            None => Err(seqalign_exception(
                DENSE_SEG_CPP,
                159,
                "CSeqalignException::eInvalidRowNumber",
                "CDense_seg::GetSeq_id(): can not get seq-id for the row requested.",
            )),
        },
        "sparse" => {
            let rows = segs.field("rows");
            let id = match row {
                0 => rows.item(0).map(|first| first.field("first-id")),
                _ => rows.item(row - 1).map(|align| align.field("second-id")),
            };
            match id {
                Some(id) => Ok(Some(id.seq_id())),
                None => Err(seqalign_exception(
                    SPARSE_SEG_CPP,
                    87,
                    "CSeqalignException::eInvalidRowNumber",
                    "CSparse_seg::GetSeq_id(): can not get seq-id for the row requested.",
                )),
            }
        }
        "dendiag" => Ok(segs
            .items()
            .find_map(|seg| seg.field("ids").item(row))
            .map(|id| id.seq_id())),
        "std" => {
            Ok(segs
                .items()
                .find_map(|seg| match seg.get("ids").and_then(|ids| ids.item(row)) {
                    Some(id) => Some(id.seq_id()),
                    None => seg.field("loc").item(row).and_then(seq_loc_id),
                }))
        }
        "disc" => {
            let at = todo.len();
            todo.extend(segs.items());
            todo[at..].reverse();
            Ok(None)
        }
        "spliced" => Ok(match row {
            0 => segs.get("product-id").map(|id| id.seq_id()),
            1 => segs.get("genomic-id").map(|id| id.seq_id()),
            _ => None,
        }),
        _ => Err(seqalign_exception(
            SEQ_ALIGN_CPP,
            391,
            "CSeqalignException::eUnsupported",
            "CSeq_align::GetSeq_id() currently does not handle this type of alignment.",
        )),
    }
}
