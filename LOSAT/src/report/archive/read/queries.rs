//! The queries of a search strategy as the FASTA text that NCBI writes and reads back: the
//! `bioseq-set` of `s_ImportQueries` written entry by entry with `CFastaOstream`, so that the
//! program reads the queries with its own FASTA reader and the command line's query options, as
//! NCBI's does (survey decision D3).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:425-438
//! ```c++
//!     } else {
//!         _ASSERT(queries.IsBioseq_set());
//!         const CBlast4_queries::TBioseq_set& bioseqs =
//!             queries.GetBioseq_set();
//!         CFastaOstream out(tmpfile->AsOutputFile(CTmpFile::eIfExists_Throw));
//!         out.SetFlag(CFastaOstream::eAssembleParts);
//!
//!         ITERATE(CBioseq_set::TSeq_set, seq_entry, bioseqs.GetSeq_set()){
//!             out.Write(**seq_entry);
//!         }
//!     }
//!
//!     const string& fname = tmpfile->GetFileName();
//!     tmpfile.Reset(new CTmpFile(fname));
//!     cmdline_args->SetInputStream(tmpfile);
//! ```
//! NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3365-3371
//! ```c++
//! void CFastaOstream::Write(const CSeq_entry& entry, const CSeq_loc* location,
//!                           bool no_scope)
//! {
//!     if (location || !no_scope) {
//!         CScope scope(*CObjectManager::GetInstance());
//!         Write(scope.AddTopLevelSeqEntry(entry), location);
//! ```
//! NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:2727-2754
//! ```c++
//! void CFastaOstream::Write(const CSeq_entry_Handle& handle,
//!                           const CSeq_loc* location)
//! {
//!     for (CBioseq_CI it(handle);  it;  ++it) {
//!         if ( !SkipBioseq(*it) ) {
//!             ...
//!             Write(*it, location);
//!         }
//!     }
//! }
//!
//! void CFastaOstream::Write(const CBioseq_Handle& handle,
//!                           const CSeq_loc* location,
//!                           const string& custom_title)
//! {
//!     WriteTitle(handle, location, custom_title);
//!     WriteSequence(handle, location);
//! }
//! ```
//! The stream's flags are the constructor's, `fInstantiateGaps | fAssembleParts | fEnableGI`, with
//! a width of 70 (sequence.cpp:2713-2720); `eAssembleParts` is a flag already set. Each entry is
//! added to a scope of its own, and its Bioseqs are written in the order of `CBioseq_CI`: depth
//! first through nested sets (probe `n_nested`, `n_nested_deep`).
//!
//! What NCBI 2.17.0 writes was read from its temporary file (`strace` of the import, probe
//! `F/f2/tools/qprobe.py` of 2026-11-20, 103 crafted query sets; the cases are the unit tests).
//! Where NCBI fails, it throws before the search and its report ends the run: the exception is
//! returned ([`RenderError::Exception`]); the file written so far is not used. What LOSAT does not
//! reproduce is refused ([`RenderError::Unsupported`]).

use std::collections::{HashMap, HashSet};

use crate::algorithm::pssm_search::annot_index::view::V;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::db_subject::bioseq_facts;
use crate::report::defline::generate_bioseq_defline;
use crate::seqid::label::{bioseq_fasta, fasta};
use crate::seqid::mask_handle::{handle_fasta_is_own, mask_handle_key};
use crate::seqid::model::ParsedSeqId;
use crate::sequence::stdaa_letters::protein_fasta_letters;
use crate::sequence::write_fasta_lines;

/// Why the queries are not written.
#[derive(Debug)]
pub enum RenderError {
    /// An exception that NCBI throws while it writes the queries; NCBI's application reports it
    /// and the run ends before any search.
    Exception(NcbiException),
    /// What LOSAT does not reproduce: the part of the queries, for a message
    /// "<part> is not supported by LOSAT's <program>".
    Unsupported(String),
}

/// `s_ImportQueries` of a `bioseq-set` (`set`, of type `Bioseq-set`): the FASTA text NCBI writes,
/// each entry of the set in turn.
pub fn render_queries(set: V<'_>) -> Result<Vec<u8>, RenderError> {
    let mut stream = FastaOstream::default();
    for entry in set.field("seq-set").items() {
        stream.write_entry(entry)?;
    }
    Ok(stream.text)
}

/// The width of the sequence lines (`CFastaOstream::SetWidth(70)` in its constructor).
const WIDTH: u32 = 70;

/// `CFastaOstream` with the constructor's flags: the text written, and `m_PreviousWholeIds`, the
/// handles of the ids written so far.
#[derive(Default)]
struct FastaOstream {
    text: Vec<u8>,
    /// The keys of the `CSeq_id_Handle`s of the ids written ([`mask_handle_key`]).
    previous_whole_ids: HashSet<Vec<u8>>,
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:2933-2944
// ```c++
//     if ( !have_range  &&  (m_Flags & fNoDupCheck) == 0) {
//         ITERATE (CBioseq::TId, id, bioseq.GetId()) {
//             CSeq_id_Handle idh = CSeq_id_Handle::GetHandle(**id);
//             pair<TSeq_id_HandleSet::iterator, bool> p
//                 = m_PreviousWholeIds.insert(idh);
//             if ( !p.second ) {
//                 NCBI_THROW(CObjmgrUtilException, eBadLocation,
//                            "Duplicate Seq-id " + (*id)->AsFastaString()
//                            + " in FASTA output");
//             }
//         }
//     }
// ```
/// The throw of a Seq-id written twice (`n_dup`: `Duplicate Seq-id lcl|Query_1 in FASTA output`).
const DUPLICATE_SEQ_ID: ThrowSite = ThrowSite {
    file: "src/objmgr/util/sequence.cpp",
    line: 2939,
    kind: "CObjmgrUtilException::eBadLocation",
    function: Some("ncbi::objects::CFastaOstream::x_WriteSeqIds"),
};

/// `CSeq_inst::GetLength()` of a Bioseq without `length` (`CSerialObject::ThrowUnassigned`; the
/// report names no function): the title's `GetInst_Length` (create_defline.cpp:732) reads it before
/// the sequence is written (`n_len_none`, `n_virtual0`).
const UNASSIGNED_LENGTH: ThrowSite = ThrowSite {
    file: "src/serial/serialobject.cpp",
    line: 306,
    kind: "CUnassignedMember::eGet",
    function: None,
};

// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:543-549
// ```c++
//         if ( cacheCoding != dataCoding || reverse ||
//              m_CaseConversion != eCaseConversion_none ) {
//             table = sx_GetConvertTable(dataCoding, cacheCoding,
//                                        reverse, m_CaseConversion);
//             if ( !table && cacheCoding != dataCoding ) {
//                 NCBI_THROW_FMT(CSeqVectorException, eCodingError,
//                                "Incompatible sequence codings: "<<
//                                dataCoding<<" -> "<<cacheCoding);
// ```
/// The throw of a Seq-data that `CSeqVector` does not convert to the letters of FASTA (`n_8na`:
/// `Incompatible sequence codings: 5 -> 1`).
const INCOMPATIBLE_CODINGS: ThrowSite = ThrowSite {
    file: "src/objmgr/seq_vector_ci.cpp",
    line: 548,
    kind: "CSeqVectorException::eCodingError",
    function: Some("ncbi::objects::CSeqVector_CI::x_FillCache"),
};

// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:58-61
// ```c++
// void ThrowOutOfRangeSeq_inst(size_t pos)
// {
//     NCBI_THROW_FMT(CSeqVectorException, eOutOfRange,
//                    "reference out of range of Seq-inst data: "<<pos);
// ```
/// The throw of a read past the end of a Seq-data (`n_len_long`: `reference out of range of
/// Seq-inst data: 6`).
const OUT_OF_RANGE: ThrowSite = ThrowSite {
    file: "src/objmgr/seq_vector_ci.cpp",
    line: 60,
    kind: "CSeqVectorException::eOutOfRange",
    function: Some("ncbi::objects::ThrowOutOfRangeSeq_inst"),
};

impl FastaOstream {
    /// `Write(entry)`: the entry's Bioseqs in the order of `CBioseq_CI`.
    fn write_entry(&mut self, entry: V<'_>) -> Result<(), RenderError> {
        let bioseqs = entry_bioseqs(entry)?;
        check_entry_ids(&bioseqs)?;
        for bioseq in bioseqs {
            self.write_bioseq(bioseq)?;
        }
        Ok(())
    }

    /// `Write(handle)`: `WriteTitle` (the ids, the title), then `WriteSequence`.
    fn write_bioseq(&mut self, bioseq: V<'_>) -> Result<(), RenderError> {
        let ids: Vec<ParsedSeqId> = bioseq
            .field("id")
            .items()
            .map(|id| id.seq_id().clone())
            .collect();
        let inst = bioseq.field("inst");
        let protein = inst.field("mol").int() == MOL_AA;
        self.write_seq_ids(bioseq, &ids, protein)?;
        let Some(length) = inst.get("length") else {
            return Err(RenderError::Exception(NcbiException::new(
                UNASSIGNED_LENGTH,
                "CSeq_inst::GetLength(): Attempt to get unassigned member NCBI-Sequence::Seq-inst.length",
            )));
        };
        let title = title_descriptor(bioseq)?;
        if title.is_none() && matches!(inst.field("repr").int(), REPR_SEG | REPR_MAP) {
            // The generator's titles of segmented and map Bioseqs (`x_SetTitleFromSegSeq`,
            // `x_SetTitleFromMap`, create_defline.cpp:4022-4028) are not modelled.
            return Err(RenderError::Unsupported(
                "a segmented or map query Bioseq without a title".into(),
            ));
        }
        self.write_seq_title(title, &ids, protein)?;
        // The reader reads a `TSeqPos` without a sign.
        let length = u32::try_from(length.int()).expect("a Seq-inst length of 32 bits");
        let letters = letters(inst, length, protein)?;
        write_fasta_lines(&letters, WIDTH, &mut self.text);
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:2946-2957
    // ```c++
    //     m_Out << '>';
    //     if (!(m_Flags & fIgnoreOriginalID) &&
    //             s_ShouldUseOriginalID(bioseq)) {
    //         string origID = s_FastaGetOriginalID(bioseq);
    //         if (! NStr::IsBlank(origID)) {
    //             m_Out << "lcl|" << origID;
    //         } else {
    //             x_WriteAsFasta(bioseq);
    //         }
    //     } else {
    //         x_WriteAsFasta(bioseq);
    //     }
    // ```
    /// `x_WriteSeqIds` without a location: the check for ids written before, `>`, then the
    /// original id of a Bioseq of local ids (`lcl|<id>`), else its ids as `x_WriteAsFasta` writes
    /// them with `fEnableGI` (`bioseq_fasta`: the GI, then the best id for FASTA, which is the
    /// same text whether or not the stream hides a prefix, as `fEnableGI` keeps them all).
    fn write_seq_ids(
        &mut self,
        bioseq: V<'_>,
        ids: &[ParsedSeqId],
        protein: bool,
    ) -> Result<(), RenderError> {
        for id in ids {
            let key = handle_key(id)?;
            if let Some(key) = key {
                if !self.previous_whole_ids.insert(key) {
                    let mut text = Vec::new();
                    fasta(id, &mut text);
                    return Err(RenderError::Exception(NcbiException::new(
                        DUPLICATE_SEQ_ID,
                        format!(
                            "Duplicate Seq-id {} in FASTA output",
                            String::from_utf8_lossy(&text)
                        ),
                    )));
                }
            }
        }
        self.text.push(b'>');
        match original_id(bioseq, ids) {
            Some(id) => {
                self.text.extend_from_slice(b"lcl|");
                self.text.extend_from_slice(id);
            }
            None => bioseq_fasta(ids, protein, &mut self.text),
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3000-3019
    // ```c++
    // void CFastaOstream::x_WriteSeqTitle(const CBioseq_Handle & bioseq_handle,
    //                                     const string& custom_title)
    // {
    //     string safe_title = (!custom_title.empty()) ? custom_title
    //         : m_Gen->GenerateDefline(bioseq_handle, x_GetTitleFlags());
    //
    //     if ( !safe_title.empty() ) {
    //         if ( !(m_Flags & fKeepGTSigns) ) {
    //             NStr::ReplaceInPlace(safe_title, ">", "_");
    //         }
    //         if (safe_title[0] != ' ') {
    //             m_Out << ' ';
    //         }
    //
    //         if ((m_Flags & fHTMLEncode) != 0) {
    //             safe_title = NStr::HtmlEncode(safe_title);
    //         }
    //         m_Out << safe_title;
    //     }
    //     m_Out << '\n';
    // }
    // ```
    /// `x_WriteSeqTitle`: `GenerateDefline` with `fFastaFormat` of the first title (a title
    /// generated where there is none), each `>` written as `_`, after a space unless it starts
    /// with one; then the line end. `fFastaFormat` changes no title of the Bioseqs written here
    /// (it drops a plasmid or WGS suffix, from descriptors that `title_descriptor` refuses), and
    /// the generator's reading past the end of a title of punctuation is approved exception 2 of
    /// PD-LOSAT-NCBI-DEFECTS (`report::defline`).
    fn write_seq_title(
        &mut self,
        title: Option<&[u8]>,
        ids: &[ParsedSeqId],
        protein: bool,
    ) -> Result<(), RenderError> {
        let generated = generate_bioseq_defline(
            title.unwrap_or_default(),
            protein,
            false,
            &bioseq_facts(ids),
        )
        // A title is a VisibleString (the reader keeps 0x20 to 0x7E), whose encoding
        // `NStr::HtmlDecode` always tells.
        .map_err(|_| {
            RenderError::Unsupported("a query title whose encoding NCBI cannot tell".into())
        })?;
        let mut safe_title = generated.text;
        for byte in &mut safe_title {
            if *byte == b'>' {
                *byte = b'_';
            }
        }
        if let Some(&first) = safe_title.first() {
            if first != b' ' {
                self.text.push(b' ');
            }
            self.text.extend_from_slice(&safe_title);
        }
        self.text.push(b'\n');
        Ok(())
    }
}

/// `Seq-inst.mol` `aa` (`CSeqVector::IsProtein`, `CDeflineGenerator`'s `m_IsAA`, `x_GetBestId`'s
/// choice of `FastaAARank`); every other molecule, `not-set` and `other` included, is written as a
/// nucleotide (probes `n_mol_*`).
const MOL_AA: i64 = 3;

/// The key of the `CSeq_id_Handle` of `id` (`None`: a handle of its own, equal to no other), or
/// the refusal of a `giim` id, whose handle depends on the ids made before it.
fn handle_key(id: &ParsedSeqId) -> Result<Option<Vec<u8>>, RenderError> {
    mask_handle_key(id).map_err(|_| RenderError::Unsupported("a giim Seq-id in the queries".into()))
}

/// The Bioseqs of an entry in the order of `CBioseq_CI` (depth first: `n_nested`,
/// `n_nested_deep`). A nested set's own descriptors are refused: what the title generator takes
/// of them is not probed past one case (`t_set_title`: a set's title is not its Bioseq's). So are
/// annotations, of a set or a Bioseq: the generator names a Bioseq without a title from the
/// features of its entry (a BioSource or Prot-ref feature, the CDS of a protein:
/// create_defline.cpp:1410, 1977, 2055), which `report::defline` does not model.
fn entry_bioseqs(entry: V<'_>) -> Result<Vec<V<'_>>, RenderError> {
    let mut out = Vec::new();
    let mut todo = vec![entry];
    while let Some(entry) = todo.pop() {
        match entry.choice() {
            ("seq", bioseq) => {
                if bioseq.get("annot").is_some() {
                    return Err(RenderError::Unsupported(
                        "a query Bioseq with annotations".into(),
                    ));
                }
                out.push(bioseq);
            }
            (_, set) => {
                if set.get("descr").is_some() || set.get("annot").is_some() {
                    return Err(RenderError::Unsupported(
                        "a nested Bioseq-set with descriptors or annotations in the queries".into(),
                    ));
                }
                let at = todo.len();
                todo.extend(set.field("seq-set").items());
                todo[at..].reverse();
            }
        }
    }
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:891-904
// ```c++
// void CTSE_Info::x_SetBioseqIds(CBioseq_Info* info)
// {
//     _ASSERT(info);
//     {{
//         CFastMutexGuard guard(m_BioseqsMutex);
//         ITERATE ( CBioseq_Info::TId, it, info->GetId() ) {
//             pair<TBioseqs::iterator, bool> ins =
//                 m_Bioseqs.insert(TBioseqs::value_type(*it, info));
//             if ( !ins.second ) {
//                 // No duplicate bioseqs in the same TSE
//                 NCBI_THROW(CObjMgrException, eAddDataError,
//                            "duplicate Bioseq id "+it->AsString()+" present in"+
//                            "\n  seq1: " + ins.first->second->IdString()+
//                            "\n  seq2: " + info->IdString());
//             }
//         }
// ```
// NCBI reference (598d8ae6): c++/src/objmgr/bioseq_info.cpp:1221-1230
// ```c++
// string CBioseq_Info::IdString(void) const
// {
//     CNcbiOstrstream os;
//     ITERATE ( TId, it, m_Id ) {
//         if ( it != m_Id.begin() )
//             os << " | ";
//         os << it->AsString();
//     }
//     return CNcbiOstrstreamToString(os);
// }
// ```
/// The throw of an id that two Bioseqs of one entry carry, or one twice (`n_dup_in_one`,
/// `d_set_dup`).
const DUPLICATE_BIOSEQ_ID: ThrowSite = ThrowSite {
    file: "src/objmgr/tse_info.cpp",
    line: 901,
    kind: "CObjMgrException::eAddDataError",
    function: Some("ncbi::objects::CTSE_Info::x_SetBioseqIds"),
};

/// `AddTopLevelSeqEntry` of the entry: its Bioseqs' ids are indexed in order by their handles, and
/// an id already indexed throws, its message spelling each id as its handle's `AsString()`. That is
/// the id's own FASTA text where [`handle_fasta_is_own`] says so (`d_set_dup`: `lcl|query_2` after
/// `lcl|Query_2`); elsewhere it can come from an id of another entry: refused.
fn check_entry_ids(bioseqs: &[V<'_>]) -> Result<(), RenderError> {
    let ids: Vec<Vec<&ParsedSeqId>> = bioseqs
        .iter()
        .map(|bioseq| bioseq.field("id").items().map(|id| id.seq_id()).collect())
        .collect();
    let mut seen: HashMap<Vec<u8>, usize> = HashMap::new();
    for (index, own) in ids.iter().enumerate() {
        for &id in own {
            let Some(key) = handle_key(id)? else {
                continue;
            };
            let Some(&first) = seen.get(&key) else {
                seen.insert(key, index);
                continue;
            };
            if !ids[first]
                .iter()
                .chain(own)
                .all(|&id| handle_fasta_is_own(id))
            {
                return Err(RenderError::Unsupported(
                    "the same Seq-id twice in one entry of the queries, spelled by another entry's id"
                        .into(),
                ));
            }
            let mut text = Vec::new();
            fasta(id, &mut text);
            return Err(RenderError::Exception(NcbiException::new(
                DUPLICATE_BIOSEQ_ID,
                format!(
                    "duplicate Bioseq id {} present in\n  seq1: {}\n  seq2: {}",
                    String::from_utf8_lossy(&text),
                    id_string(&ids[first]),
                    id_string(own)
                ),
            )));
        }
    }
    Ok(())
}

/// `CBioseq_Info::IdString` of ids whose handles spell them as they are.
fn id_string(ids: &[&ParsedSeqId]) -> String {
    let mut text = Vec::new();
    for (n, id) in ids.iter().enumerate() {
        if n > 0 {
            text.extend_from_slice(b" | ");
        }
        fasta(id, &mut text);
    }
    String::from_utf8_lossy(&text).into_owned()
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:2757-2810
// ```c++
// static string s_FastaGetOriginalID (const CBioseq& seq)
//
// {
//     FOR_EACH_SEQDESC_ON_BIOSEQ (it, seq) {
//         const CSeqdesc& desc = **it;
//         if (! desc.IsUser()) continue;
//         if (! desc.GetUser().IsSetType()) continue;
//         const CUser_object& usr = desc.GetUser();
//         const CObject_id& oi = usr.GetType();
//         if (! oi.IsStr()) continue;
//         const string& type = oi.GetStr();
//         if (! NStr::EqualNocase(type, "OrginalID") && ! NStr::EqualNocase(type, "OriginalID")) continue;
//         FOR_EACH_USERFIELD_ON_USEROBJECT (uitr, usr) {
//             const CUser_field& fld = **uitr;
//             if (FIELD_IS_SET_AND_IS(fld, Label, Str)) {
//                 const string &label_str = GET_FIELD(fld.GetLabel(), Str);
//                 if (! NStr::EqualNocase(label_str, "LocalId")) continue;
//                 if (fld.IsSetData() && fld.GetData().IsStr()) {
//                     return fld.GetData().GetStr();
//                 }
//             }
//         }
//     }
//
//     return "";
// }
//
// static bool s_ShouldUseOriginalID (const CBioseq& seq)
// {
//     FOR_EACH_SEQID_ON_BIOSEQ (id_itr, seq) {
//         const CSeq_id& sid = **id_itr;
//         switch (sid.Which()) {
//             case CSeq_id::e_Local:
//                 break;
//             case CSeq_id::e_General:
//                 {
//                     const CDbtag& dbtag = sid.GetGeneral();
//                     if (dbtag.IsSetDb()) {
//                         const string& db = dbtag.GetDb();
//                         if (! NStr::EqualNocase(db, "TMSMART") &&
//                             ! NStr::EqualNocase(db, "BankIt") &&
//                             ! NStr::EqualNocase(db, "NCBIFILE")) {
//                             return false;
//                         }
//                     }
//                 }
//                 break;
//             default:
//                 return false;
//         }
//     }
//
//     return true;
// }
// ```
/// The original id that `x_WriteSeqIds` writes as `lcl|<id>`: that of the first `OriginalID` (or
/// `OrginalID`) user object with a `LocalId` string field, for a Bioseq whose ids are all local or
/// of the general databases TMSMART, BankIt and NCBIFILE (probes `n_orig_*`); none when it is
/// blank (`NStr::IsBlank`: white space only).
fn original_id<'a>(bioseq: V<'a>, ids: &[ParsedSeqId]) -> Option<&'a [u8]> {
    let should_use = ids.iter().all(|id| match id {
        ParsedSeqId::Local(_) => true,
        ParsedSeqId::General(dbtag) => [&b"TMSMART"[..], b"BankIt", b"NCBIFILE"]
            .iter()
            .any(|db| dbtag.db.eq_ignore_ascii_case(db)),
        _ => false,
    });
    if !should_use {
        return None;
    }
    let found = descriptors(bioseq)
        .filter_map(|(variant, user)| (variant == "user").then_some(user))
        .filter(|user| {
            matches!(user.field("type").choice(), ("str", kind)
                if is_original_id_type(kind.str()))
        })
        .find_map(|user| {
            user.field("data").items().find_map(|field| {
                match (field.field("label").choice(), field.field("data").choice()) {
                    (("str", label), ("str", text))
                        if label.str().eq_ignore_ascii_case(b"LocalId") =>
                    {
                        Some(text.str())
                    }
                    _ => None,
                }
            })
        })?;
    // `NStr::IsBlank`: `isspace` of every byte, in the C locale.
    (!found
        .iter()
        .all(|&byte| matches!(byte, b' ' | b'\t' | b'\n' | 0x0B | 0x0C | b'\r')))
    .then_some(found)
}

/// Whether a user object's type is that of an original id (`OrginalID`, `OriginalID`, any case).
fn is_original_id_type(kind: &[u8]) -> bool {
    kind.eq_ignore_ascii_case(b"OriginalID") || kind.eq_ignore_ascii_case(b"OrginalID")
}

/// The Bioseq's descriptors in order: their variant and value.
fn descriptors<'a>(bioseq: V<'a>) -> impl Iterator<Item = (&'static str, V<'a>)> {
    bioseq
        .get("descr")
        .into_iter()
        .flat_map(|descr| descr.items())
        .map(|desc| desc.choice())
}

/// The first `title` descriptor, which the title generator takes (`n_title_two`), after the
/// refusal of the descriptors that would change the generated title in ways `report::defline`
/// does not model: the generator reads `molinfo` (tech prefixes and suffixes: `n_molinfo_tsa`
/// gives `TSA: tsa title, mRNA sequence`), `source`, `genbank` and `embl` (keywords), `comment`,
/// `pdb`, and the `Unverified`, `Unreviewed` and `AutodefOptions` user objects
/// (create_defline.cpp:840-1135). A user object whose type is a string other than `CFastaReader`
/// (the reader's own, which NCBI's exports carry) or an original id is refused with them; one
/// typed by a number is read by neither. The other descriptors are read by neither.
fn title_descriptor(bioseq: V<'_>) -> Result<Option<&[u8]>, RenderError> {
    let mut title = None;
    for (variant, desc) in descriptors(bioseq) {
        match variant {
            "title" => {
                title.get_or_insert(desc.str());
            }
            "user" => {
                if let ("str", kind) = desc.field("type").choice() {
                    let kind = kind.str();
                    if kind != b"CFastaReader" && !is_original_id_type(kind) {
                        return Err(RenderError::Unsupported(format!(
                            "a query's user object of type {}",
                            String::from_utf8_lossy(kind)
                        )));
                    }
                }
            }
            "molinfo" | "source" | "genbank" | "embl" | "comment" | "pdb" => {
                return Err(RenderError::Unsupported(format!(
                    "a query's {variant} descriptor"
                )));
            }
            _ => {}
        }
    }
    Ok(title)
}

/// A segment of the sequence map: a Seq-data and its length, or a gap.
enum Piece<'a> {
    Data(V<'a>, u32),
    Gap(u32),
}

impl Piece<'_> {
    fn length(&self) -> u32 {
        match *self {
            Piece::Data(_, length) | Piece::Gap(length) => length,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:708-714
// ```c++
//     if ( !m_Seg && pos == x_GetSize() ) {
//         // it's ok to position to the very end
//         return;
//     }
//     if ( !m_Seg || pos<m_Seg.GetPosition() || pos>=m_Seg.GetEndPosition() ) {
//         NCBI_THROW_FMT(CSeqVectorException, eDataError,
//                        "CSeqVector_CI: cannot locate segment at "<<pos);
// ```
/// The throw of a map whose segments end before the Bioseq's length (`dl_sum_long`: `cannot
/// locate segment at 9`). NCBI's kind is `eSeqDataError`, the text of `eDataError`.
const NO_SEGMENT: ThrowSite = ThrowSite {
    file: "src/objmgr/seq_vector_ci.cpp",
    line: 713,
    kind: "CSeqVectorException::eSeqDataError",
    function: Some("ncbi::objects::CSeqVector_CI::x_UpdateSeg"),
};

/// `WriteSequence` without a location: the letters of the Bioseq's `CSeqVector` in
/// `eCoding_Iupac`, `Ncbieaa` for a protein (sequence.cpp:3322-3361), read through
/// `CSeqVector_CI` up to the Bioseq's length, which is the vector's size (`GetBioseqLength`); gaps
/// instantiated as `N` (`X` for a protein: `n_raw_nodata`, `p_virtual`, `p_delta_gap`). A segment
/// that passes the length is read up to it (`dl_sum_short`); segments that end before it throw
/// where they end.
fn letters(inst: V<'_>, length: u32, protein: bool) -> Result<Vec<u8>, RenderError> {
    let pieces = sequence_map(inst, length)?;
    let target = if protein { NCBIEAA } else { IUPACNA };
    let mut out = Vec::new();
    let mut position = 0u32;
    for piece in pieces {
        if position >= length {
            break;
        }
        let Some(end) = position.checked_add(piece.length()) else {
            return Err(RenderError::Unsupported(
                "a delta query Bioseq whose segments pass 2^32 letters".into(),
            ));
        };
        let needed = end.min(length) - position;
        match piece {
            Piece::Gap(_) => out.resize(
                out.len() + needed as usize,
                if protein { b'X' } else { b'N' },
            ),
            Piece::Data(data, piece_length) => {
                data_letters(data, piece_length, needed, target, &mut out)?
            }
        }
        position = end;
    }
    if position < length {
        return Err(RenderError::Exception(NcbiException::new(
            NO_SEGMENT,
            format!("CSeqVector_CI: cannot locate segment at {position}"),
        )));
    }
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:219-291
// ```c++
//     if ( inst.IsSetSeq_data() ) {
//         if ( !inst.GetSeq_data().IsGap() ) {
//             x_Add(inst.GetSeq_data(), inst.GetLength());
//         }
//         else {
//             // split Seq-data
//             x_AddGap(inst.GetLength(), false, inst.GetSeq_data());
//         }
//     }
//     else if ( inst.IsSetExt() ) {
//         const CSeq_ext& ext = inst.GetExt();
//         switch (ext.Which()) {
//         case CSeq_ext::e_Seg:
//             ITERATE ( CSeq_ext::TSeg::Tdata, iter, ext.GetSeg().Get() ) {
//                 x_Add(**iter);
//             }
//             break;
//         case CSeq_ext::e_Ref:
//             x_Add(ext.GetRef());
//             break;
//         case CSeq_ext::e_Delta:
//             ITERATE ( CSeq_ext::TDelta::Tdata, iter, ext.GetDelta().Get() ) {
//                 x_Add(**iter);
//             }
//             break;
//         case CSeq_ext::e_Map:
//             if ( m_SeqLength == kInvalidSeqPos ) {
//                 NCBI_THROW(CSeqMapException, eDataError,
//                        "CSeq_inst.ext.map with unset CSeq_inst.length");
//             }
//             x_AddGap(m_SeqLength, false);
//             break;
//         default:
//             //### Not implemented
//             NCBI_THROW(CSeqMapException, eUnimplemented,
//                        "CSeq_ext::??? -- not implemented");
//         }
//     }
//     else if ( inst.GetRepr() == CSeq_inst::eRepr_virtual ) {
//         // Virtual sequence -- no data, no segments
//         // The total sequence is gap
//         if ( m_SeqLength == kInvalidSeqPos ) {
//             m_SeqLength = 0;
//         }
//         x_AddGap(m_SeqLength, false);
//     }
//     else if ( inst.GetRepr() == CSeq_inst::eRepr_map ) {
//         // Virtual map sequence -- no data, no segments
//         // The total sequence is gap
//         if ( m_SeqLength == kInvalidSeqPos ) {
//             NCBI_THROW(CSeqMapException, eDataError,
//                        "CSeq_inst.repr map with unset CSeq_inst.length");
//         }
//         x_AddGap(m_SeqLength, false);
//     }
//     else if ( inst.GetRepr() != CSeq_inst::eRepr_not_set &&
//               inst.IsSetLength() && inst.GetLength() != 0 ) {
//         // split seq-data
//         x_AddGap(inst.GetLength(), false);
//     }
//     else {
//         if ( inst.GetRepr() != CSeq_inst::eRepr_not_set ) {
//             NCBI_THROW(CSeqMapException, eDataError,
//                        "CSeq_inst.repr of sequence without data "
//                        "should be not_set");
//         }
//         if ( inst.IsSetLength() && inst.GetLength() != 0 ) {
//             NCBI_THROW(CSeqMapException, eDataError,
//                        "CSeq_inst.length of sequence without data "
//                        "should be 0");
//         }
//         x_AddGap(0, false);
//     }
// ```
/// The Bioseq's sequence map (`CSeqMap(const CSeq_inst&)`) for a Bioseq whose length is set: its
/// Seq-data over the length, whatever its representation (`r_*`); else its extension (a delta's
/// literals in order, a literal without Seq-data or with `gap` data being a gap; a map's gap over
/// the length); else a gap over the length (`virtual`, `map`, or another representation of a
/// length other than 0: `n_raw_nodata`); else nothing, or the map's throw. Segments of other
/// Bioseqs (`seg`, `ref`, a delta's locations) are refused: NCBI fails to find those Bioseqs
/// (`n_seg`, `n_delta_loc`: `Cannot resolve lcl|Other: unknown`), in a message LOSAT does not
/// reproduce.
fn sequence_map(inst: V<'_>, length: u32) -> Result<Vec<Piece<'_>>, RenderError> {
    if let Some(data) = inst.get("seq-data") {
        return Ok(vec![data_piece(data, length)]);
    }
    if let Some(ext) = inst.get("ext") {
        return match ext.choice() {
            ("delta", delta) => delta
                .items()
                .map(|item| match item.choice() {
                    ("literal", literal) => {
                        let length = seq_length(literal.field("length").int())?;
                        Ok(match literal.get("seq-data") {
                            Some(data) => data_piece(data, length),
                            None => Piece::Gap(length),
                        })
                    }
                    _ => Err(far_pointer()),
                })
                .collect(),
            ("map", _) => Ok(vec![Piece::Gap(length)]),
            _ => Err(far_pointer()),
        };
    }
    match inst.field("repr").int() {
        REPR_VIRTUAL | REPR_MAP => Ok(vec![Piece::Gap(length)]),
        REPR_NOT_SET if length == 0 => Ok(vec![Piece::Gap(0)]),
        REPR_NOT_SET => Err(seq_map_error(
            "CSeq_inst.length of sequence without data should be 0",
        )),
        _ if length != 0 => Ok(vec![Piece::Gap(length)]),
        _ => Err(seq_map_error(
            "CSeq_inst.repr of sequence without data should be not_set",
        )),
    }
}

/// The throw of `CSeqMap(const CSeq_inst&)` for a Bioseq without data whose representation or
/// length says it has some (`m_raw_len0_nodata`, `m_notset_len4_nodata`).
fn seq_map_error(message: &str) -> RenderError {
    let line = if message.ends_with("not_set") {
        281
    } else {
        286
    };
    RenderError::Exception(NcbiException::new(
        ThrowSite {
            file: "src/objmgr/seq_map.cpp",
            line,
            kind: "CSeqMapException::eDataError",
            function: Some("ncbi::objects::CSeqMap::CSeqMap"),
        },
        message,
    ))
}

/// The refusal of a segment of another Bioseq.
fn far_pointer() -> RenderError {
    RenderError::Unsupported(
        "a query Bioseq made of other Bioseqs (seg, ref or a delta's location)".into(),
    )
}

/// `Seq-inst.repr` values.
const REPR_NOT_SET: i64 = 0;
const REPR_VIRTUAL: i64 = 1;
const REPR_SEG: i64 = 3;
const REPR_MAP: i64 = 7;

/// A literal's length (`TSeqPos`).
fn seq_length(length: i64) -> Result<u32, RenderError> {
    u32::try_from(length)
        .map_err(|_| RenderError::Unsupported("a delta literal length outside TSeqPos".into()))
}

/// A Seq-data of `length` letters: a gap when it is `gap` data.
fn data_piece(data: V<'_>, length: u32) -> Piece<'_> {
    match data.choice() {
        ("gap", _) => Piece::Gap(length),
        _ => Piece::Data(data, length),
    }
}

/// `CSeq_data::E_Choice` of the codings: the numbers of NCBI's messages.
const IUPACNA: u32 = 1;
const IUPACAA: u32 = 2;
const NCBI2NA: u32 = 3;
const NCBI4NA: u32 = 4;
const NCBIEAA: u32 = 8;
const NCBISTDAA: u32 = 10;

/// The number of a Seq-data variant (`CSeq_data::E_Choice`).
fn coding(variant: &str) -> u32 {
    match variant {
        "iupacna" => IUPACNA,
        "iupacaa" => IUPACAA,
        "ncbi2na" => NCBI2NA,
        "ncbi4na" => NCBI4NA,
        "ncbi8na" => 5,
        "ncbipna" => 6,
        "ncbi8aa" => 7,
        "ncbieaa" => NCBIEAA,
        "ncbipaa" => 9,
        "ncbistdaa" => NCBISTDAA,
        other => unreachable!("a Seq-data variant of the table, not {other}"),
    }
}

/// NCBI4na to IUPACna as `CSeqVector` converts (code 0, a gap, is `N`: `n_4na_codes`).
const NA4_TO_IUPACNA: &[u8; 16] = b"NACMGRSVTWYHKDBN";
/// NCBI2na to IUPACna, high bits first.
const NA2_TO_IUPACNA: &[u8; 4] = b"ACGT";

// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:484-494
// ```c++
// void CSeqVector_CI::x_UpdateCacheUp(TSeqPos pos)
// {
//     _ASSERT(pos < x_GetSize());
//
//     TSeqPos segEnd = m_Seg.GetEndPosition();
//     _ASSERT(pos >= m_Seg.GetPosition() && pos < segEnd);
//
//     TSeqPos cache_size = min(kCacheSize, segEnd - pos);
//     x_FillCache(pos, cache_size);
// ```
/// The first `needed` letters of a segment of `length` letters of a Seq-data, filled as
/// `CSeqVector_CI` fills its cache: from the start of the segment, 1024 letters at a time
/// (`kCacheSize`) but not past the segment's end, while a fill starts before `needed` (the
/// iterator's end: a fill of a segment that passes the Bioseq's length can pass it). The fill converts the
/// coding or throws (`Incompatible sequence codings`: a coding of the other molecule, NCBI8na,
/// NCBIpna, NCBI8aa, NCBIpaa; probes `n_8na`, `n_pna`, `p_8aa`, `p_paa`, `n_mol_*`, `p_mol_*`); each
/// fill checks its end against the data (`copy_8bit_any`, `copy_4bit_any`, `copy_2bit_any` of
/// include/objmgr/impl/seq_vector_cvt.hpp:47-105: the end `e` of the fill throws when `e`, `e / 2`
/// or `e / 4` passes the data's bytes). A packed coding's fill that ends inside the byte past the
/// data reads memory the input does not decide: refused.
fn data_letters(
    data: V<'_>,
    length: u32,
    needed: u32,
    target: u32,
    out: &mut Vec<u8>,
) -> Result<(), RenderError> {
    if needed == 0 {
        return Ok(());
    }
    let (variant, value) = data.choice();
    let source = coding(variant);
    let convertible = matches!(
        (source, target),
        (IUPACNA | NCBI2NA | NCBI4NA, IUPACNA) | (NCBIEAA | IUPACAA | NCBISTDAA, NCBIEAA)
    );
    if !convertible {
        return Err(RenderError::Exception(NcbiException::new(
            INCOMPATIBLE_CODINGS,
            format!("Incompatible sequence codings: {source} -> {target}"),
        )));
    }
    let bytes = match variant {
        "iupacna" | "iupacaa" | "ncbieaa" => value.str(),
        _ => value.octets(),
    };
    let per_byte: u64 = match source {
        NCBI2NA => 4,
        NCBI4NA => 2,
        _ => 1,
    };
    let held = bytes.len() as u64 * per_byte;
    // The fill that holds the first letter past the data, if it is made: it throws, or reads past
    // the data's end.
    let start = held / KCACHE_SIZE * KCACHE_SIZE;
    if u64::from(length) > held && start < u64::from(needed) {
        let end = (start + KCACHE_SIZE).min(u64::from(length));
        if end / per_byte > bytes.len() as u64 {
            return Err(RenderError::Exception(NcbiException::new(
                OUT_OF_RANGE,
                format!("reference out of range of Seq-inst data: {end}"),
            )));
        }
        return Err(RenderError::Unsupported(
            "a packed query Seq-data that ends inside its last byte before the Bioseq's length (NCBI reads past the data)"
                .into(),
        ));
    }
    let length = needed as usize;
    match source {
        IUPACNA | IUPACAA | NCBIEAA => out.extend_from_slice(&bytes[..length]),
        NCBISTDAA => out.extend(protein_fasta_letters(&bytes[..length])),
        NCBI4NA => out.extend((0..length).map(|i| {
            let byte = bytes[i / 2];
            NA4_TO_IUPACNA[usize::from(if i % 2 == 0 { byte >> 4 } else { byte & 15 })]
        })),
        NCBI2NA => out.extend((0..length).map(|i| {
            let shift = 6 - 2 * (i % 4);
            NA2_TO_IUPACNA[usize::from((bytes[i / 4] >> shift) & 3)]
        })),
        _ => unreachable!("a convertible coding"),
    }
    Ok(())
}

/// `kCacheSize` of `CSeqVector_CI` (seq_vector_ci.cpp:56).
const KCACHE_SIZE: u64 = 1024;

#[cfg(test)]
#[path = "queries_tests.rs"]
mod tests;
