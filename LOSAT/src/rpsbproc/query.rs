//! The queries of one archive with their domain hits: `CNcbiCdAlignProcessor::ProcessBlastResults`
//! for protein queries (each query `Bioseq` with the alignments that follow it in the list), the
//! query's label, title and length (`FillSequenceFromBioseq`), and the masks of the results
//! (`ConvertBlastMaskListToSeqLocInfoVector`, which rpsbproc reads but does not print).

use super::archive::{Align, Mask, MaskLocation, Query, Segs, SeqEntry};
use super::data::{compute_pseudo_pssm_id, CddData};
use super::domains::{calculate, parse_align, AlignIndex, DomAlign};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::seqid::defline_set::{compare, Sic};
use crate::seqid::model::{ObjectId, ParsedSeqId};

/// A failure while the queries of an archive are processed.
#[derive(Debug)]
pub(crate) enum ProcessError {
    /// An exception of NCBI's serial objects thrown while a session is processed
    /// (`CInvalidChoiceSelection::eFail` from a getter of a CHOICE that holds another variant,
    /// `CUnassignedMember::eGet` from a getter of an unset member; rpsbproc: two `Error:` lines,
    /// exit 3).
    Exception(NcbiException),
    /// What LOSAT does not reproduce (an explicit refusal).
    Unsupported(String),
}

/// The `CInvalidChoiceSelection` message of `C<type>::Get<variant>()` on a CHOICE of
/// `module::type` that holds `current`.
// NCBI reference (598d8ae6): c++/src/serial/exception.cpp:122-150
// ```c++
// CInvalidChoiceSelection::CInvalidChoiceSelection(
//     const CDiagCompileInfo& diag_info,
//     const CSerialObject* object,
//     size_t currentIndex, size_t mustBeIndex,
//     const char* const names[], size_t namesCount,
//     EDiagSev severity)
//         : CSerialException(diag_info, 0, "", severity, 0)
// {
//     CNcbiOstrstream msg;
//     const CChoiceTypeInfo* type = 0;
//     if ( object ) {
//         type = dynamic_cast<const CChoiceTypeInfo*>(object->GetThisTypeInfo());
//     }
//     const char* cur_name = GetName(currentIndex, names, namesCount);
//     const char* req_name = GetName(mustBeIndex, names, namesCount);
//     if ( type ) {
//         msg << "C" << SPrintIdentifier(type->GetAccessName())
//             << "::Get" << SPrintIdentifier(req_name) << "()";
//         msg << ": Invalid choice selection: "
//             << type->GetAccessModuleName() <<"::"<< type->GetAccessName()
//             << '.' << cur_name;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/serial/exception.cpp:86-118 (`SPrintIdentifier`: `-` is
// `_`, the first letter and each one after `.` upper-cased)
// (NCBI's rpsbproc printed `CSeq_id::GetGeneral(): Invalid choice selection:
// NCBI-Seqloc::Seq-id.local` for a profile database with `lcl|` ids, PS16a probe 2.)
// (The exception is thrown by the generated `C<type>_Base::ThrowInvalidSelection`, whose place
// LOSAT does not have: `ThrowSite::generated`.)
pub(crate) fn choice_selection(
    module: &str,
    type_name: &str,
    required: &str,
    current: &str,
) -> ProcessError {
    ProcessError::Exception(NcbiException::new(
        ThrowSite::generated("CInvalidChoiceSelection::eFail"),
        format!(
            "C{}::Get{}(): Invalid choice selection: {module}::{type_name}.{current}",
            identifier(type_name),
            identifier(required)
        ),
    ))
}

/// `SPrintIdentifier` of a type or member name without `.`: `-` read as `_`, the first letter
/// upper-cased.
fn identifier(name: &str) -> String {
    let mut out = String::new();
    let mut capitalize = true;
    for c in name.chars() {
        let c = if c == '-' { '_' } else { c };
        if capitalize {
            out.extend(c.to_uppercase());
            capitalize = false;
        } else {
            out.push(c);
        }
    }
    out
}

/// The `CUnassignedMember` exception of `C<type>::Get<member>()` on a SEQUENCE of
/// `module::type` whose `member` is not set.
// NCBI reference (598d8ae6): c++/src/serial/serialobject.cpp:266-306
// ```c++
// void CSerialObject::ThrowUnassigned(TMemberIndex index,
//                                     const char* file_name,
//                                     int file_line) const
// {
//     if (x_GetVerifyData() == eSerialVerifyData_Yes) {
//         ...
//         CNcbiOstrstream s;
//         if ( true ) {
//             // make class name
//             s << "C" << SPrintIdentifier(type->GetAccessName());
//         }
//         if ( !member_name.empty() ) {
//             // make method name
//             s << "::Get" << SPrintIdentifier(member_name) << "()";
//         }
//         s << ": Attempt to get unassigned member "
//           << type->GetAccessModuleName() <<"::"<< type->GetAccessName() << '.';
//         if ( !member_name.empty() ) {
//             s << member_name;
//         ...
//         NCBI_THROW(CUnassignedMember,eGet,CNcbiOstrstreamToString(s));
// ```
// (`x_GetVerifyData` is `eSerialVerifyData_Yes` unless `SERIAL_VERIFY_DATA_GET` says otherwise,
// which LOSAT's rpsbproc rejects: `reader_variable`.)
pub(crate) fn unassigned_member(module: &str, type_name: &str, member: &str) -> ProcessError {
    ProcessError::Exception(NcbiException::new(
        UNASSIGNED_MEMBER,
        format!(
            "C{}::Get{}(): Attempt to get unassigned member {module}::{type_name}.{member}",
            identifier(type_name),
            identifier(member)
        ),
    ))
}

// NCBI reference (598d8ae6): c++/src/serial/serialobject.cpp:317-320
// ```c++
// void CSerialObject::ThrowUnassigned(TMemberIndex index) const
// {
//     ThrowUnassigned(index, 0, 0);
// }
// ```
// NCBI reference (598d8ae6): c++/src/serial/datatool/classstr.cpp:849-855
// ```c++
//                     code.Methods(inl) <<
//                         "    if (!CanGet"<< mem_cName<<"()) {\n"
//                         "        ThrowUnassigned("<<member_index;
// #if 0
//                     code.Methods(inl) << ", __FILE__, __LINE__";
// #endif
//                     code.Methods(inl) << ");\n"
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:805-810
// ```c++
// const char* find_match(char        lsep,
//                        char        rsep,
//                        const char* start,
//                        const char* stop)
// {
//     if (*(stop - 1) != rsep) return stop;
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:836-842
// ```c++
//     const char* end_str = find_match('(', ')',
//                                      m_CurrFunctName,
//                                      m_CurrFunctName + len);
//     if (end_str == m_CurrFunctName + len) {
//         // Missing '('
//         return;
//     }
// ```
// The generated getters call `ThrowUnassigned(index)`, so the place is the `NCBI_THROW` of
// serialobject.cpp:306 (`file_name` and `file_line` are 0). Its function is a `const` method,
// whose `__PRETTY_FUNCTION__` ends with `const`, not `)`: `ParseCurrFunctName` names no
// function then.
const UNASSIGNED_MEMBER: ThrowSite = ThrowSite {
    file: "src/serial/serialobject.cpp",
    line: 306,
    kind: "CUnassignedMember::eGet",
    function: None,
};

/// A protein query of the output (`TDomQuery`) with its hits.
#[derive(Debug)]
pub(crate) struct DomQuery {
    /// `m_strAccession`: the label of the query's best id.
    pub(crate) accession: Vec<u8>,
    /// `m_strTitle`.
    pub(crate) title: Vec<u8>,
    /// `m_uiSeqLen`.
    pub(crate) length: u32,
    pub(crate) aligns: Vec<DomAlign>,
    pub(crate) index: AlignIndex,
}

/// `ConvertBlastMaskListToSeqLocInfoVector` as far as it can fail: the first location of each
/// mask that has locations must be a `packed-int`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:264-300
// ```c++
// void ConvertBlastMaskListToSeqLocInfoVector(const TBlast4MaskList &bmlist, TSeqLocInfoVector &dst)
// {
//
// 	TBlast4MaskList::const_iterator iterBM = bmlist.begin(), iterBMEnd = bmlist.end();
//
// 	const CSeq_id *pLastId = nullptr;
// 	TSeqLocInfoVector::iterator iterDstMaskedRegions = dst.end();
// 	while (iterBMEnd != iterBM)
// 	{
// 		const CBlast4_mask::TLocations &locs = (*iterBM)->GetLocations();	//list< CRef< CSeq_loc > >
//
// 		if (!locs.empty())	//should always be true, expect only one
// 		{
// 			int rf = (*iterBM)->GetFrame();
// 			if (rf > 0)
// 				rf = READINGFRAME::Idx2Id(rf - 1);	//CSeqLocInfo::ETranslationFrame uses 0, 1, 2, 3, -1, -2, -3
// 			const CSeq_loc& loc = **locs.begin();
// 			const CSeq_id *thisId = loc.GetId();	//should not be null
//
// 			if (nullptr == pLastId || !thisId->Match(*pLastId))
// 			{
//
// 				// -- create a new TMaskedQueryRegions
// 				iterDstMaskedRegions = dst.emplace(dst.end());
// 				pLastId = thisId;
// 			}
//
// 			const list< CRef< CSeq_interval > > & ints = loc.GetPacked_int().Get();
//
// 			for (auto v : ints)
// 				iterDstMaskedRegions->emplace_back(new CSeqLocInfo(v.GetNCPointer(), rf));
// 		}
//
// 		++iterBM;
//
// 	}
// }
// ```
// The masks are kept with the queries (`m_vecMaskedRegions`), which the text output never prints;
// what is observable is the exception of a location that is not a `packed-int`. A location
// without one id (`GetId` gives no id: `null`, `feat`, `mix`, `equiv`, `bond`, or intervals of
// different ids) is dereferenced by NCBI: refused.
pub(crate) fn check_masks(masks: &[Mask]) -> Result<(), ProcessError> {
    for mask in masks {
        let Some(first) = &mask.first else {
            continue;
        };
        match first {
            MaskLocation::PackedInt(ids) => {
                if ids.is_empty() || ids.iter().any(|id| id != &ids[0]) {
                    return Err(ProcessError::Unsupported(
                        "a mask whose packed-int location is empty or names more than one Seq-id is not supported by LOSAT's rpsbproc".to_string(),
                    ));
                }
            }
            MaskLocation::Other(name @ ("int" | "whole" | "empty" | "pnt" | "packed-pnt")) => {
                return Err(choice_selection(
                    "NCBI-Seqloc",
                    "Seq-loc",
                    "packed-int",
                    name,
                ));
            }
            MaskLocation::Other(name) => {
                return Err(ProcessError::Unsupported(format!(
                    "a mask location of the Seq-loc variant {name} is not supported by LOSAT's rpsbproc"
                )));
            }
        }
    }
    Ok(())
}

/// `CSeq_id::GetStringDescr(bioseq, eFormat_BestWithVersion)` and the clean-up of
/// `FillSequenceFromBioseq`: the label of the best id, the text after its first `:` and then
/// after its first `|`.
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2292-2345
// ```c++
// string CSeq_id::GetStringDescr(const CBioseq& bioseq, EStringFormat fmt)
// {
//     ...
//     bool is_na            = bioseq.GetInst().GetMol() != CSeq_inst::eMol_aa;
//     CRef<CSeq_id> best_id = FindBestChoice(bioseq.GetId(),
//                                            is_na ? CSeq_id::FastaNARank
//                                            : CSeq_id::FastaAARank);
//     switch (fmt) {
//     ...
//     case eFormat_BestWithVersion:
//         // eBestWithVersion produces only the 'best' accession name, with
//         // its version indicator
//         if (best_id.NotEmpty()) {
//             string label;
//             best_id->GetLabel(&label, eDefault, fLabel_Version);
//             return label;
//         }
//         break;
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:1030-1046
// ```c++
// 	if (TSequence::e_IsLocal == dst.m_iValid)	//replace accession with label
// 	{
// 		// -- check original fasta defline from reader
// 		dst.m_OriDefline.clear();
// 		ObtainUserDeflineFromBioseq(bioseq, dst.m_OriDefline);
//
// 		if (dst.m_strAccession.empty() && !localId.empty())
// 			dst.m_strAccession = move(localId);
// 	}
//
// 	// -- sanitize accession
// 	size_t dlm_pos = dst.m_strAccession.find(':');
// 	if (string::npos != dlm_pos)
// 		dst.m_strAccession = dst.m_strAccession.substr(dlm_pos + 1);
// 	dlm_pos = dst.m_strAccession.find('|');
// 	if (string::npos != dlm_pos)
// 		dst.m_strAccession = dst.m_strAccession.substr(dlm_pos + 1);
// }
// ```
// (The best id's label is never empty when the Bioseq has an id, so the local id's label is
// never used.)
fn query_accession(query: &Query) -> Vec<u8> {
    let is_na = query.mol != 3;
    let score = if is_na {
        crate::seqid::rank::fasta_na_score
    } else {
        crate::seqid::rank::fasta_aa_score
    };
    let mut label = Vec::new();
    if let Some(best) = crate::seqid::rank::find_best(&query.ids, score) {
        crate::seqid::label::both_version(best, &mut label);
    }
    if let Some(at) = label.iter().position(|&b| b == b':') {
        label.drain(..=at);
    }
    if let Some(at) = label.iter().position(|&b| b == b'|') {
        label.drain(..=at);
    }
    label
}

/// `FillSequenceFromBioseq` for the members that rpsbproc prints: the label, the title and the
/// length of a protein query.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:952-983
// ```c++
// void FillSequenceFromBioseq(const CBioseq &bioseq, TSequence &dst)
// {
// 	if (bioseq.IsSetInst() && bioseq.CanGetInst())
// 	{
// 		dst.m_iValid = TSequence::e_Current;	//by default
// 		const CSeq_inst & seqInst = bioseq.GetInst();
// 		dst.m_strSeqData = Get1LetterSeqData(seqInst);
//
//
// 		if (seqInst.IsSetMol() && seqInst.CanGetMol())
// 			dst.m_bIsNa = (CSeq_inst::eMol_aa != seqInst.GetMol());
// 			//dst.m_iMolType = seqInst.GetMol();
//
// 		if (seqInst.IsSetLength() && seqInst.CanGetLength())
// 			dst.m_uiSeqLen = seqInst.GetLength();
// 		else
// 			dst.m_uiSeqLen = dst.m_strSeqData.size();
// 	}
// 	...
// 	dst.m_strTitle = GetDescrTitleFromBioseq(bioseq);
// 	if (0 == dst.m_iTaxId) dst.m_iTaxId = bioseq.GetTaxId();
//
// 	const CBioseq::TId &ids = bioseq.GetId();
// 	string localId;
//     CBioseq::TId pack_ids;  //list to collect ids to be packed
//
// 	dst.m_strAccession = CSeq_id::GetStringDescr(bioseq, CSeq_id::eFormat_BestWithVersion);
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:248-253
// ```c++
// 	// -- somehow the CSeqportUtil::Convert may append an extra residue to the result (bug?), so check
// 	SeqLen_t ori_len = rSeqInst.GetLength(), end_len = dst.size();
// 	if (ori_len < end_len)
// 	{
// 		dst.erase(ori_len);
// 	}
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:297-317
// ```c++
// const string& GetDescrTitleFromBioseq(const CBioseq& rBioseq)
// {
// 	if (rBioseq.IsSetDescr() && rBioseq.CanGetDescr())
// 	{
// 		...
// 				CSeqdesc::E_Choice c = rSeqDesc.Which();
// 				if (CSeqdesc::e_Title == c)
// 				{
// 					return rSeqDesc.GetTitle();
// 				}
// 		...
// 	if (rBioseq.IsSetAnnot() && rBioseq.CanGetAnnot())
// 	{
// ```
// `Get1LetterSeqData` converts the residues, then asks for `inst.length` whether it is set or
// not: `GetLength` of an unset member throws `CUnassignedMember`, which ends rpsbproc as NCBI's
// does (exit 3, after the `SESSION` line). The residues serve the site output of item PS16b. A
// query without a title descriptor but with `annot` takes its title from a feature (refused),
// and a `CFastaReader` DefLine that is not a string is refused (`Query::unread_descriptor`).
// `GetTaxId` (biodata_core.cpp:977) reads the `org` and `source` descriptors but never throws,
// and the taxid is not printed, so those descriptors are only read past.
fn fill_sequence(query: &Query) -> Result<(Vec<u8>, Vec<u8>, u32), ProcessError> {
    // `SeqDataTo1LtrSeq` (biodata_core.cpp:149-197) converts the protein encodings with
    // `CSeqportUtil::Convert`; LOSAT reads the one-letter ones and NCBIstdaa's 28 codes and
    // refuses the others (before the length, as the conversion comes first).
    match &query.seq_data {
        None => {}
        Some(("ncbieaa" | "iupacaa", _)) => {}
        Some(("ncbistdaa", codes)) if codes.iter().all(|&code| code < 28) => {}
        Some((name, _)) => {
            return Err(ProcessError::Unsupported(format!(
                "a protein query whose Seq-data is {name} is not supported by LOSAT's rpsbproc"
            )));
        }
    }
    let Some(length) = query.length else {
        return Err(unassigned_member("NCBI-Sequence", "Seq-inst", "length"));
    };
    if let Some(what) = query.unread_descriptor {
        let local = query
            .ids
            .iter()
            .any(|id| matches!(id, ParsedSeqId::Local(_)));
        // `ObtainUserDeflineFromBioseq` runs for a query with a local id only.
        if local {
            return Err(ProcessError::Unsupported(format!(
                "a query Bioseq with {what} is not supported by LOSAT's rpsbproc"
            )));
        }
    }
    let title = match &query.title {
        Some(title) => title.clone(),
        None if query.has_annot => {
            return Err(ProcessError::Unsupported(
                "a query Bioseq without a title descriptor but with annot (NCBI takes the title from a feature) is not supported by LOSAT's rpsbproc".to_string(),
            ));
        }
        None => Vec::new(),
    };
    Ok((query_accession(query), title, length))
}

/// `x_GetPssmId`: the PSSM-ID of a subject id (`ProcessCDAlign` keeps the hits whose PSSM-ID
/// is not 0).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_biodata.cpp:175-191
// ```c++
// PssmId_t CNcbiCdAlignProcessor::x_GetPssmId(const CSeq_id &seqid, string & acxn) const
// {
//
// 	const CSeq_id::TGeneral& rDbTag = seqid.GetGeneral();  //CDbtag
//
// 	const CDbtag::TTag& rObjectId = rDbTag.GetTag();  //CObject_id
//
// 	acxn.clear();
// 	if (rObjectId.Which() == CObject_id::e_Id)	//correct pssmid
// 		return rObjectId.GetId();
// 	else if (rObjectId.Which() == CObject_id::e_Str)
// 	{
// 		acxn = rObjectId.GetStr();
// 		return TDomSrcCount::ComputePseudoPssmId(acxn);
// 	}
// 	return 0;
// }
// ```
// (The `int` id is returned as `PssmId_t`, an `unsigned int`.)
fn pssm_id(id: &ParsedSeqId) -> Result<u32, ProcessError> {
    let ParsedSeqId::General(dbtag) = id else {
        return Err(choice_selection(
            "NCBI-Seqloc",
            "Seq-id",
            "general",
            id.alternative_name(),
        ));
    };
    Ok(match &dbtag.tag {
        ObjectId::Id(value) => *value as u32,
        ObjectId::Str(text) => compute_pseudo_pssm_id(text),
    })
}

/// `ProcessCDQuery` of a protein query: its hits (`ProcessCDAlign`) and their classification.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_biodata.cpp:195-232
// ```c++
// void CNcbiCdAlignProcessor::ProcessCDAlign(const list<CRef<CSeq_align> > &rAligns, vector<TDomSeqAlignment> &dst) const
// {
// 	for (list<CRef<CSeq_align> > :: const_iterator iterSrc = rAligns.begin(); iterSrc != rAligns.end(); ++iterSrc)
// 	{
//
// 		TDomSeqAlignment newalign;
//
// 		CRef<CSeq_id> hit_id = ParseAlignSegs((*iterSrc)->GetSegs(), newalign);
//
// 		if (!hit_id.IsNull())
// 		{
// 			string acxn;
// 			PssmId_t uiPssmId = x_GetPssmId(*hit_id, acxn);
// 			if (uiPssmId > 0)	//pseudo or not
// 			{
//
// 				newalign.m_uiPssmId = uiPssmId;
//
// 				ParseAlignScores(**iterSrc, newalign);
//
// 				dst.emplace_back(move(newalign));
// 			}
// 		}
// 	}
// }
//
// void CNcbiCdAlignProcessor::ProcessCDQuery(const list<ncbi::CRef<ncbi::objects::CSeq_align> > &rAligns, TDomQuery &dst, vector<PssmId_t> *missed) const
// {
//
// 	ProcessCDAlign(rAligns, dst.m_vecAlignments);
// 	...
// 	if (!dst.m_bIsNa)
// 	{
// 		Calculate(dst.m_vecAlignments, dst.m_dimSplitAligns[0], missed_pssmids);
// 	}
// ```
fn process_cd_query(
    aligns: &[Align],
    query: &mut DomQuery,
    data: &CddData,
) -> Result<(), ProcessError> {
    for align in aligns {
        let Segs::Denseg(seg) = &align.segs else {
            unreachable!("the caller refuses other segments");
        };
        // (`process_blast_results` refuses a Dense-seg without ids before this.)
        if seg.ids.len() < 2 || seg.starts.len() < 2 * seg.lens.len() {
            return Err(ProcessError::Unsupported(
                "a Dense-seg with fewer than 2 ids or fewer starts than 2 per segment (NCBI reads past them) is not supported by LOSAT's rpsbproc".to_string(),
            ));
        }
        let pssm_id = pssm_id(&seg.ids[1])?;
        if pssm_id == 0 {
            continue;
        }
        if let Some(name) = align.score_of_other_type {
            return Err(ProcessError::Unsupported(format!(
                "a Seq-align whose {name} score has a value of another type is not supported by LOSAT's rpsbproc"
            )));
        }
        // A NaN E-value is not ordered with the others: `TDomAlignOrderByEValue` is then no strict
        // weak order, and NCBI's `std::sort` of the hits is undefined once it partitions (more
        // than 16 hits; with 16 or fewer its insertion sort gives a defined order). The archives
        // of rpsblast hold none (`NOT-A-NUMBER` is only a REAL's text form); LOSAT refuses every
        // NaN E-value rather than keep two rules by the number of hits, and a NaN bit score,
        // whose comparisons and printing it does not port.
        if align.e_value.is_some_and(f64::is_nan) || align.bit_score.is_some_and(f64::is_nan) {
            return Err(ProcessError::Unsupported(
                "a Seq-align whose e_value or bit_score is not a number is not supported by LOSAT's rpsbproc".to_string(),
            ));
        }
        query.aligns.push(parse_align(seg, align, pssm_id));
    }
    query.index = calculate(&mut query.aligns, data);
    Ok(())
}

/// `ProcessBlastResults` of protein queries: each query with the alignments that follow in the
/// list while their first id matches the query's local id (else its first id).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_biodata.cpp:308-366
// ```c++
// void CNcbiCdAlignProcessor::ProcessBlastResults(list<TDomQuery> &dst, const list<CRef<CSeq_entry> > &qseqs, const TSeqLocInfoVector & masks, const list< CRef<CSeq_align > > &aligns, int gcode, vector<PssmId_t> *missed) const
// {
//
// 	if (!qseqs.empty())
// 	{
// 		list< CRef<CSeq_align> > :: const_iterator iterAlign = aligns.begin(), iterAlignEnd = aligns.end();
// 		...
// 		for (list<CRef<CSeq_entry> > :: const_iterator iterSeq = qseqs.begin(), iterSeqEnd = qseqs.end(); iterSeqEnd != iterSeq; ++iterSeq)
// 		{
// 			const CBioseq &seq = (*iterSeq)->GetSeq();	//offline sequence, assume to be seq, no more check
// 			list<TDomQuery> :: iterator iterDQ = dst.emplace(dst.end());
//
// 			FillSequenceFromBioseq(seq, *iterDQ);
//
// 			const CSeq_id* curr_seq_id = seq.GetLocalId();
// 			if (nullptr == curr_seq_id)	//added because Bioseq from command line rpsblast may not contain a local seq-id.
// 				curr_seq_id = seq.GetFirstId();
// 			...
// 			list< CRef< CSeq_align > > lstAligns;
//
// 			while (iterAlignEnd != iterAlign)
// 			{
// 				const CSeq_id &aligned = (*iterAlign)->GetSeq_id(0);
// 				if (!curr_seq_id->Match(aligned))
// 				{
// 					if (!lstAligns.empty())
// 						ProcessCDQuery(lstAligns, *iterDQ, missed);
//
// 					goto labelNextSeq;	//to next sequence
// 				}
//
// 				// -- collect aligns for current sequence
// 				lstAligns.emplace_back(*iterAlign);
//
// 				++iterAlign;
// 			}
//
// 			if (!lstAligns.empty())
// 				ProcessCDQuery(lstAligns, *iterDQ, missed);
//
// 		labelNextSeq:;
// 		}
// ```
// (`CSeq_id::Match` is `Compare(...) == e_YES`, `seqid::defline_set::compare`. The masks of
// each query are not printed. A Seq-align other than a `Dense-seg` comes from a nucleotide
// search, item PS16c: refused.)
pub(crate) fn process_blast_results(
    entries: &[SeqEntry],
    aligns: &[Align],
    data: &CddData,
) -> Result<Vec<DomQuery>, ProcessError> {
    let mut queries = Vec::with_capacity(entries.len());
    let mut next = 0usize;
    for entry in entries {
        let SeqEntry::Seq(query) = entry else {
            return Err(choice_selection("NCBI-Seqset", "Seq-entry", "seq", "set"));
        };
        let (accession, title, length) = fill_sequence(query)?;
        let mut dom_query = DomQuery {
            accession,
            title,
            length,
            aligns: Vec::new(),
            index: AlignIndex::default(),
        };
        // NCBI reference (598d8ae6): c++/src/objects/seq/Bioseq.cpp:271-279,327-343
        // ```c++
        // const CSeq_id* CBioseq::GetFirstId() const
        // {
        //     // If no ids for Bioseq, return 0 -- should not happen
        //     if (GetId().empty()) {
        //         return 0;
        //     }
        //
        //     return *GetId().begin();
        // }
        // ...
        // static int s_BestLocalRank(const CRef<CSeq_id>& id)
        // {
        //     if (id.Empty() || !id->IsLocal()) {
        //         return kMax_Int;
        //     }
        //     return id->BestRankScore();
        // }
        //
        // const CSeq_id* CBioseq::GetLocalId() const
        // {
        //     CRef<CSeq_id> id = FindBestChoice(GetId(), &s_BestLocalRank);
        //     if (id.NotEmpty() && id->IsLocal()) {
        //         return &*id;
        //     }
        //     return NULL;
        // }
        // ```
        let current = crate::seqid::rank::find_best(&query.ids, |id| {
            if matches!(id, ParsedSeqId::Local(_)) {
                crate::seqid::rank::best_rank_score(id)
            } else {
                i32::MAX
            }
        })
        .or(query.ids.first());
        let start = next;
        while next < aligns.len() {
            let Segs::Denseg(seg) = &aligns[next].segs else {
                let Segs::Other(name) = &aligns[next].segs else {
                    unreachable!()
                };
                return Err(ProcessError::Unsupported(format!(
                    "a Seq-align with {name} segments (a nucleotide search's, item PS16c) is not supported by LOSAT's rpsbproc"
                )));
            };
            let Some(aligned) = seg.ids.first() else {
                return Err(ProcessError::Unsupported(
                    "a Dense-seg without ids (NCBI's CSeq_align::GetSeq_id throws) is not supported by LOSAT's rpsbproc".to_string(),
                ));
            };
            let Some(current) = current else {
                return Err(ProcessError::Unsupported(
                    "a query Bioseq without ids (NCBI dereferences its missing id) is not supported by LOSAT's rpsbproc".to_string(),
                ));
            };
            if compare(current, aligned) != Sic::Yes {
                break;
            }
            next += 1;
        }
        if next > start {
            process_cd_query(&aligns[start..next], &mut dom_query, data)?;
        }
        queries.push(dom_query);
    }
    Ok(queries)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqid::model::{Dbtag, TextseqId};

    // serial/exception.cpp:86-150 (quoted at `choice_selection`); NCBI's rpsbproc printed the first
    // one for a profile database with `lcl|` ids (PS16a probe 2).
    #[test]
    fn choice_selection_messages() {
        let message = |error| match error {
            ProcessError::Exception(exception) => exception.get_msg().to_string(),
            ProcessError::Unsupported(message) => panic!("{message}"),
        };
        assert_eq!(
            message(choice_selection(
                "NCBI-Seqloc",
                "Seq-id",
                "general",
                "local"
            )),
            "CSeq_id::GetGeneral(): Invalid choice selection: NCBI-Seqloc::Seq-id.local"
        );
        assert_eq!(
            message(choice_selection("NCBI-Blast4", "Blast4-cutoff", "e-value", "raw-score")),
            "CBlast4_cutoff::GetE_value(): Invalid choice selection: NCBI-Blast4::Blast4-cutoff.raw-score"
        );
        assert_eq!(
            message(choice_selection(
                "NCBI-Seqloc",
                "Seq-loc",
                "packed-int",
                "int"
            )),
            "CSeq_loc::GetPacked_int(): Invalid choice selection: NCBI-Seqloc::Seq-loc.int"
        );
    }

    // serialobject.cpp:266-306 (quoted at `unassigned_member`); NCBI's rpsbproc printed it for a
    // query without `inst.length` (PS16a review oracle O3-9).
    #[test]
    fn unassigned_member_message() {
        let ProcessError::Exception(exception) =
            unassigned_member("NCBI-Sequence", "Seq-inst", "length")
        else {
            panic!("an exception");
        };
        assert_eq!(exception.site().kind, "CUnassignedMember::eGet");
        assert_eq!(
            exception.get_msg(),
            "CSeq_inst::GetLength(): Attempt to get unassigned member NCBI-Sequence::Seq-inst.length"
        );
    }

    // Seq_id.cpp:2292-2326 and biodata_core.cpp:1030-1046 (quoted at `query_accession`): the
    // label of the best id, after its first `:`, then after its first `|`.
    #[test]
    fn query_labels() {
        let query = |ids: Vec<ParsedSeqId>| Query {
            ids,
            mol: 3,
            ..Query::default()
        };
        let local = ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()));
        assert_eq!(query_accession(&query(vec![local.clone()])), b"Query_1");
        let general = ParsedSeqId::General(Dbtag {
            db: b"CDD".to_vec(),
            tag: ObjectId::Id(5),
        });
        assert_eq!(query_accession(&query(vec![general])), b"5");
        let refseq = ParsedSeqId::Other(TextseqId {
            accession: Some(b"NP_414783".to_vec()),
            version: Some(1),
            ..TextseqId::default()
        });
        // FastaAARank prefers the RefSeq id to the local one.
        assert_eq!(query_accession(&query(vec![local, refseq])), b"NP_414783.1");
        assert_eq!(query_accession(&query(vec![ParsedSeqId::Gi(123)])), b"123");
        assert_eq!(query_accession(&query(Vec::new())), b"");
    }
}
