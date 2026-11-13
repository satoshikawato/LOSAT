//! The queries of one archive with their domain hits: `CNcbiCdAlignProcessor::ProcessBlastResults`
//! (each query `Bioseq` with the alignments that follow it in the list; a nucleotide query's
//! hits classified per reading frame), the query's label, title, length and residues
//! (`FillSequenceFromBioseq`), and the masks of the results
//! (`ConvertBlastMaskListToSeqLocInfoVector`, which rpsbproc reads but does not print).

use super::archive::{Align, Mask, MaskLocation, Query, QueryRowId, Segs, SeqEntry};
use super::data::{compute_pseudo_pssm_id, CddData};
use super::domains::{
    calculate, parse_align, parse_std_align, sort_reading_frames, AlignIndex, DomAlign,
};
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
    /// A `std::exception` that leaves `Run` (`CSimpleException` of `THROW_SIMPLE`): its
    /// `what()`, which NCBI's application layer reports (exit 3).
    Std(String),
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

/// A query of the output (`TDomQuery`) with its hits.
#[derive(Debug)]
pub(crate) struct DomQuery {
    /// `m_strAccession`: the label of the query's best id.
    pub(crate) accession: Vec<u8>,
    /// `m_strTitle`.
    pub(crate) title: Vec<u8>,
    /// `m_uiSeqLen`.
    pub(crate) length: u32,
    /// `m_bIsNa`: `inst.mol` is not `aa`.
    pub(crate) is_na: bool,
    /// `m_strSeqData`: the residues (or bases) as one-letter codes (`Get1LetterSeqData`); empty
    /// when the query has neither `seq-data` nor literals of `ext delta` with data.
    pub(crate) residues: Vec<u8>,
    pub(crate) aligns: Vec<DomAlign>,
    /// `m_dimSplitAligns`: the classification of the hits of each reading frame (index 0 only
    /// for a protein query).
    pub(crate) frames: [AlignIndex; 6],
    /// `_TOfflDomQuery::m_dimTranslated` of a query of the XML input: the residues that its
    /// alignments give each reading frame (`-` elsewhere), which the sites of a nucleotide query
    /// read in place of `TranslateAll` of its bases; `None` for a query of an archive.
    pub(crate) translated: Option<Box<[Vec<u8>; 6]>>,
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

/// `FillSequenceFromBioseq` for the members that rpsbproc prints: the label, the title, the
/// length and the residues (or bases) of a query.
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
// does (exit 3, after the `SESSION` line). The residues serve the site output (PS16b; a
// nucleotide query's are translated, PS16c). A
// query without a title descriptor but with `annot` takes its title from a feature (refused),
// and a `CFastaReader` DefLine that is not a string is refused (`Query::unread_descriptor`).
// `GetTaxId` (biodata_core.cpp:977) reads the `org` and `source` descriptors but never throws,
// and the taxid is not printed, so those descriptors are only read past.
fn fill_sequence(query: &Query) -> Result<(Vec<u8>, Vec<u8>, u32, Vec<u8>), ProcessError> {
    // `SeqDataTo1LtrSeq` (biodata_core.cpp:149-197) converts the encodings with
    // `CSeqportUtil::Convert`, which LOSAT ports for every encoding but NCBI8aa and NCBIstdaa
    // codes from 28 (refused, before the length, as the conversion comes first; NCBIpaa and
    // NCBIpna throw NCBI's exception), for the `seq-data` and for each literal of `ext delta`
    // (`Get1LetterSeqData` converts each, in order).
    let residues = one_letter_residues(query).map_err(|(what, name)| {
        if name == UNCONVERTIBLE {
            return ProcessError::Std(UNCONVERTIBLE.to_string());
        }
        ProcessError::Unsupported(format!(
            "a query whose {what} is {name} is not supported by LOSAT's rpsbproc"
        ))
    })?;
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
    let mut residues = residues;
    // `Get1LetterSeqData` cuts the letters to `inst.length` (`GetLength`).
    residues.truncate(length as usize);
    Ok((query_accession(query), title, length, residues))
}

/// The kind of a `Seq-data` after `SeqDataTo1LtrSeq`: 0 for the protein encodings, 1 for the
/// nucleotide ones.
type LetterType = u8;

/// `SeqDataTo1LtrSeq`: the letters and their kind, or the name of an encoding that LOSAT does
/// not convert. (`gap` data, type -1, is refused as a `seq-data`, as in PS16a: the reader does
/// not check its grammar; the archive reader leaves the literals of `gap` data out.)
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:149-197
// ```c++
// int SeqDataTo1LtrSeq(const CSeq_data &rSeqData, string &dst)
// {
// 	dst.clear();
// 	int iDataType = rSeqData.Which();
//
// 	switch (iDataType)
// 	{
// 	case CSeq_data::e_Ncbieaa:
//
// 		dst = rSeqData.GetNcbieaa().Get();
//
// 		iDataType = 0;
// 		break;
// 	case CSeq_data::e_Ncbi8aa:
// 	case CSeq_data::e_Ncbipaa:
// 	case CSeq_data::e_Iupacaa:
// 	case CSeq_data::e_Ncbistdaa:	//all other aa format
// 		{
// 			CSeq_data objEaaData;
// 			CSeqportUtil::Convert(rSeqData, &objEaaData, CSeq_data::e_Ncbieaa);
// 			dst = objEaaData.GetNcbieaa().Get();
// 		}
// 		iDataType = 0;
// 		break;
// 	case CSeq_data::e_Iupacna:
// 		dst = rSeqData.GetIupacna().Get();
//
//
// 		iDataType = 1;
// 		break;
// 	case CSeq_data::e_Ncbi2na:
// 	case CSeq_data::e_Ncbi4na:
// 	case CSeq_data::e_Ncbi8na:
// 	case CSeq_data::e_Ncbipna:
// 		{
// 			CSeq_data objIupacData;
// 			CSeqportUtil::Convert(rSeqData, &objIupacData, CSeq_data::e_Iupacna);
//
//
// 			dst = objIupacData.GetIupacna().Get();
//
// 		}
// 		iDataType = 1;
// 		break;
// 	default:
// 		iDataType = -1;
// 	}
// 	return iDataType;
// }
// ```
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:367-374,651-658,777-784
// ```c++
// SIZE_TYPE CSeqConvert_imp::x_Convert2naToIupacna
// ...
//     return convert_1_to_4(src, pos, length, dst, C2naToIupacna::GetTable());
// ...
// SIZE_TYPE CSeqConvert_imp::x_Convert4naToIupacna
// ...
//     return convert_1_to_2(src, pos, length, dst, C4naToIupacna::GetTable());
// ...
// SIZE_TYPE CSeqConvert_imp::x_Convert8naToIupacna
// ...
//     return convert_1_to_1(src, pos, length, dst, C8naToIupacna::GetTable());
// ```
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1504-1507
// ```c++
// const Uint1 C8naToIupacna::scm_Table[256] = {
//     0x4e, 0x41, 0x43, 0x4d, 0x47, 0x52, 0x53, 0x56,
//     0x54, 0x57, 0x59, 0x48, 0x4b, 0x44, 0x42, 0x4e,
//     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
// ```
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:886-893,950-958
// ```c++
// SIZE_TYPE CSeqConvert_imp::x_ConvertIupacaaToEaa(...)
// {
//     return Subseq(src, CSeqUtil::e_Iupacaa, pos, length, dst);
// }
// ...
//     return convert_1_to_1(src, pos, length, dst, CStdaaToEaa::GetTable());
// ```
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1822-1826
// ```c++
// const Uint1 CStdaaToEaa::scm_Table[256] = {
//     0x2d, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47,
//     0x48, 0x49, 0x4b, 0x4c, 0x4d, 0x4e, 0x50, 0x51,
//     0x52, 0x53, 0x54, 0x56, 0x57, 0x58, 0x59, 0x5a,
//     0x55, 0x2a, 0x4f, 0x4a, 0xff, 0xff, 0xff, 0xff,
// ```
// The tables `C2naToIupacna` and `C4naToIupacna` (sequtil_tables.cpp:498, 1094) spell the
// bases of each byte from its high bits: NCBI2na codes 0 to 3 are `ACGT`, NCBI4na codes 0 to
// 15 are `NACMGRSVTWYHKDBN` (code 0, a gap, is `N`); NCBI8na codes 16 and up are the byte
// 0xff. A whole number of bytes is converted (`AdjustLength` keeps every base of the data), and
// `Get1LetterSeqData` cuts the result to `inst.length`. iupacaa to ncbieaa copies the letters;
// ncbistdaa codes 0 to 27 are the characters `-ABCDEFGHIKLMNPQRSTVWXYZU*OJ` and the other
// codes, 0xff, are refused as in PS16a, and so is NCBI8aa (PS16a). NCBIpaa and NCBIpna have no
// conversion in `CSeqportUtil`: NCBI's `runtime_error` (`UNCONVERTIBLE`).
fn one_letter(data: &(&'static str, Vec<u8>)) -> Result<(LetterType, Vec<u8>), &'static str> {
    const STDAA_TO_EAA: &[u8; 28] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
    const NA2: &[u8; 4] = b"ACGT";
    const NA4: &[u8; 16] = b"NACMGRSVTWYHKDBN";
    let (name, bytes) = data;
    Ok(match *name {
        "ncbieaa" | "iupacaa" => (0, bytes.clone()),
        "ncbistdaa" => {
            if bytes.iter().any(|&code| code >= 28) {
                return Err(name);
            }
            (
                0,
                bytes
                    .iter()
                    .map(|&code| STDAA_TO_EAA[code as usize])
                    .collect(),
            )
        }
        "iupacna" => (1, bytes.clone()),
        "ncbi2na" => (
            1,
            bytes
                .iter()
                .flat_map(|&b| [6u8, 4, 2, 0].map(|shift| NA2[usize::from((b >> shift) & 3)]))
                .collect(),
        ),
        "ncbi4na" => (
            1,
            bytes
                .iter()
                .flat_map(|&b| [NA4[usize::from(b >> 4)], NA4[usize::from(b & 15)]])
                .collect(),
        ),
        "ncbi8na" => (
            1,
            bytes
                .iter()
                .map(|&b| if b < 16 { NA4[usize::from(b)] } else { 0xff })
                .collect(),
        ),
        "ncbipna" | "ncbipaa" => return Err(UNCONVERTIBLE),
        _ => return Err(name),
    })
}

/// `CSeqportUtil`'s exception for an encoding without a conversion (NCBIpna, NCBIpaa): a
/// `std::runtime_error`, reported by `CNcbiApplication` (oracle `pna` of PS16c: `Error: (106.7)
/// Application's execution failed: Requested conversion not implemented`, exit 3).
// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:2029-2035
// ```c++
//     if(to_code == CSeq_data::e_not_set || from_code == CSeq_data::e_not_set) {
//         throw std::runtime_error("to_code or from_code not set");
//     }
//     if ( s_SeqDataToSeqUtil[to_code]  == CSeqUtil::e_not_set  ||
//          s_SeqDataToSeqUtil[from_code] == CSeqUtil::e_not_set ) {
//         throw runtime_error("Requested conversion not implemented");
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:1947-1959
// ```c++
// static CSeqUtil::TCoding s_SeqDataToSeqUtil[] = {
//     CSeqUtil::e_not_set,
//     CSeqUtil::e_Iupacna,
//     CSeqUtil::e_Iupacaa,
//     CSeqUtil::e_Ncbi2na,
//     CSeqUtil::e_Ncbi4na,
//     CSeqUtil::e_Ncbi8na,
//     CSeqUtil::e_not_set,
//     CSeqUtil::e_Ncbi8aa,
//     CSeqUtil::e_Ncbieaa,
//     CSeqUtil::e_not_set,
//     CSeqUtil::e_Ncbistdaa
// };
// ```
const UNCONVERTIBLE: &str = "Requested conversion not implemented";

/// `Get1LetterSeqData` before its cut to `inst.length`: the `seq-data`, else the literals of
/// `ext delta` of the kind of the first one with letters, one after another; or what to refuse
/// (the place and the encoding).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:199-246
// ```c++
// string Get1LetterSeqData(const CSeq_inst& rSeqInst)
// {
// 	string dst;
// 	if (rSeqInst.IsSetSeq_data() && rSeqInst.CanGetSeq_data())
// 	{
// 		SeqDataTo1LtrSeq(rSeqInst.GetSeq_data(), dst);
// 	}
// 	else if (rSeqInst.IsSetExt() && rSeqInst.CanGetExt())
// 	{
// 		const CSeq_ext &rSeqExt = rSeqInst.GetExt();
//
// 		switch (rSeqExt.Which())
// 		{
// 		case CSeq_ext::e_not_set:
// 		case CSeq_ext::e_Seg:
// 		case CSeq_ext::e_Ref:
// 		case CSeq_ext::e_Map:
// 			break;
// 		case CSeq_ext::e_Delta:
// 			{
// 				const CDelta_ext::Tdata &rDeltaExt_data = rSeqExt.GetDelta().Get();	//typedef list< CRef< CDelta_seq > > Tdata;
// 				int litType = -1;
// 				for (CDelta_ext::Tdata::const_iterator iter = rDeltaExt_data.begin(), iterEnd = rDeltaExt_data.end(); iterEnd != iter; ++iter)
// 				{
// 					switch ((*iter)->Which())
// 					{
// 					case CDelta_seq::e_Literal:
// 						{
// 							const CSeq_literal &rSeqLit = (*iter)->GetLiteral();
// 							if (rSeqLit.IsSetSeq_data() && rSeqLit.CanGetSeq_data())
// 							{
// 								string tmp(k_strEmptyString);
// 								int stype = SeqDataTo1LtrSeq(rSeqLit.GetSeq_data(), tmp);
// 								if (stype >= 0)	//valid
// 								{
// 									if (litType < 0) litType = stype;
//
// 									if (litType == stype)	//same aa or na
// 										dst.append(tmp);
// 								}
// 							}
// 						}
// 						break;
// 					}
// 				}
// 			}
// 		}
// 	}
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
fn one_letter_residues(query: &Query) -> Result<Vec<u8>, (&'static str, &'static str)> {
    if let Some(data) = &query.seq_data {
        return Ok(one_letter(data).map_err(|name| ("Seq-data", name))?.1);
    }
    let mut residues = Vec::new();
    let mut literal_type = None;
    for data in &query.delta_literals {
        let (kind, letters) =
            one_letter(data).map_err(|name| ("Seq-data of a delta literal", name))?;
        if *literal_type.get_or_insert(kind) == kind {
            residues.extend_from_slice(&letters);
        }
    }
    Ok(residues)
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

/// `ProcessCDQuery`: the query's hits (`ProcessCDAlign`) and their classification, of all of
/// them for a protein query, of each reading frame's for a nucleotide query.
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
// 	else
// 	{
// 		vector<size_t> rfIndice[READINGFRAME::TOTAL_RFS];
// 		CCdAlignProcessor::SortReadingFrames(dst.m_vecAlignments, rfIndice, dst.m_uiSeqLen);
// 		for (int i = 0; i < READINGFRAME::TOTAL_RFS; ++i)
// 		{
//
// 			if (!rfIndice[i].empty())
// 			{
// 				vecMissed.clear();
//
// 				Calculate(dst.m_vecAlignments, rfIndice[i], dst.m_dimSplitAligns[i], missed_pssmids);
//
// 			}
// 		}
// 	}
// }
// ```
// (`ParseAlignSegs` is `domains::parse_align` and `domains::parse_std_align`; it takes the
// hit's id from the second id of a Dense-seg or of the first Std-seg, and reads past `ids` when
// there are fewer: refused. The reading frames of a nucleotide query's hits are the classes of
// `domains::sort_reading_frames`; its exception for a query of length 0 is NCBI's, after the
// hits are read.)
fn process_cd_query(
    aligns: &[Align],
    query: &mut DomQuery,
    data: &CddData,
) -> Result<(), ProcessError> {
    for align in aligns {
        let (dom, hit) = match &align.segs {
            Segs::Denseg(seg) => {
                // (`process_blast_results` refuses a Dense-seg without ids before this.)
                if seg.ids.len() < 2 || seg.starts.len() < 2 * seg.lens.len() {
                    return Err(ProcessError::Unsupported(
                        "a Dense-seg with fewer than 2 ids or fewer starts than 2 per segment (NCBI reads past them) is not supported by LOSAT's rpsbproc".to_string(),
                    ));
                }
                (None, &seg.ids[1])
            }
            Segs::Std(std) => {
                // (`process_blast_results` finds a first id in the segments before this, so
                // there is at least one.)
                let Some(hit) = &std.hit_id else {
                    return Err(ProcessError::Unsupported(
                        "a Seq-align whose first Std-seg has fewer than 2 ids (NCBI reads past them) is not supported by LOSAT's rpsbproc".to_string(),
                    ));
                };
                let dom =
                    parse_std_align(&std.segs, align, 0).map_err(ProcessError::Unsupported)?;
                (Some(dom), hit)
            }
            Segs::Other(_) => unreachable!("the caller refuses other segments"),
        };
        let pssm_id = pssm_id(hit)?;
        if pssm_id == 0 {
            continue;
        }
        if let Some(name) = align.score_of_other_type {
            return Err(ProcessError::Unsupported(format!(
                "a Seq-align whose {name} score has a value of another type is not supported by LOSAT's rpsbproc"
            )));
        }
        // (A NaN E-value or bit score, `NOT-A-NUMBER` in the text, is kept: `domains::calculate`
        // sorts the hits as NCBI's `std::sort` does and refuses only where that sort would read
        // outside its vector; the DOMAINS line prints it `-nan`, `asn_text_reader::NCBI_NAN`.)
        let dom = match (dom, &align.segs) {
            (Some(mut dom), _) => {
                dom.pssm_id = pssm_id;
                dom
            }
            (None, Segs::Denseg(seg)) => parse_align(seg, align, pssm_id),
            (None, _) => unreachable!("a Dense-seg"),
        };
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1809,1836
        // ```c++
        // 	(*m_ostream) << pType << DELIMIT << (rec.pAlign->m_uiPssmId >= INVALIDPSSMID ? 0 : rec.pAlign->m_uiPssmId) << DELIMIT << rec.pAlign->m_iFrom + COORDSBASE << DELIMIT << rec.pAlign->m_iTo + COORDSBASE << DELIMIT << ...
        // 	(*m_ostream) << HITTYPE_CLUSTER << DELIMIT << (rec.pCdInfo->m_uiClusterPssmId >= INVALIDPSSMID ? 0 : rec.pCdInfo->m_uiClusterPssmId) << DELIMIT << rec.pAlign->m_iFrom + COORDSBASE << DELIMIT << rec.pAlign->m_iTo + COORDSBASE << DELIMIT << ...
        // ```
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:205-206
        // ```c++
        // 			CSegSet thisSegs;
        // 			thisSegs.AddSeg(iter->pAlign->m_iFrom, iter->pAlign->m_iTo);
        // ```
        // A kept hit of a nucleotide query whose range starts or ends at 2^31 - 1: the `int`
        // sums of its DOMAINS and SUPERFAMILIES lines overflow (undefined; PS16c review oracle
        // `plusmax`: both 0.5.1 builds print -2147483648), and so may the overlap sets of
        // `Calculate` that start from the range (`GetTotalResidues`' `to - from + 1`, the
        // `+ 1` of `Cross` and `Merge`): refused, whatever the strand (PD-LOSAT-NCBI-DEFECTS
        // 1.16). A protein query's hit keeps PS16a's widened sums (`print_dom_line`).
        if query.is_na && (dom.from == i32::MAX || dom.to == i32::MAX) {
            return Err(ProcessError::Unsupported(
                "a hit of a nucleotide query that starts or ends at position 2147483647 (NCBI's int overflows) is not supported by LOSAT's rpsbproc".to_string(),
            ));
        }
        query.aligns.push(dom);
    }
    if !query.is_na {
        let all: Vec<usize> = (0..query.aligns.len()).collect();
        query.frames[0] =
            calculate(&mut query.aligns, &all, data).map_err(ProcessError::Unsupported)?;
        return Ok(());
    }
    let frames = sort_reading_frames(&query.aligns, query.length)
        .map_err(|what| ProcessError::Std(what.to_string()))?;
    for (frame, selected) in frames.iter().enumerate() {
        if !selected.is_empty() {
            query.frames[frame] =
                calculate(&mut query.aligns, selected, data).map_err(ProcessError::Unsupported)?;
        }
    }
    Ok(())
}

/// `CSeq_align::GetSeq_id(0)`: the id of the query row of an alignment, or what LOSAT refuses
/// (NCBI's exceptions when no segment names one, and the location variants whose id LOSAT does
/// not look up).
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:317-358
// ```c++
// const CSeq_id& CSeq_align::GetSeq_id(TDim row) const
// {
//     switch (GetSegs().Which()) {
//     case C_Segs::e_Denseg:
//         return GetSegs().GetDenseg().GetSeq_id(row);
//     ...
//     case C_Segs::e_Std:
//         {{
//             // If segments have different number of rows, this will try
//             // to find the segment with enough rows to get the id.
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
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:390-398
// ```c++
//     default:
//         NCBI_THROW(CSeqalignException, eUnsupported,
//                    "CSeq_align::GetSeq_id() currently does not handle "
//                    "this type of alignment.");
//     }
//
//     NCBI_THROW(CSeqalignException, eInvalidRowNumber,
//                "CSeq_align::GetSeq_id(): "
//                "can not get seq-id for the row requested.");
// ```
// (A `null` location has no id; `empty`, `whole`, `int` and `pnt` have theirs. The exceptions
// end rpsbproc as its other `CException`s of a session do: PS16c oracles `noids2` and
// `packed`, the two `Error:` lines, exit 3. `GetSeq_id` is a `const` method, so the report names
// no function, as `UNASSIGNED_MEMBER`'s. The `packed` segments are read for their grammar only:
// NCBI's `ParseAlignSegs` has a branch for them, which this exception keeps it from reaching.
// `dendiag`, `disc`, `spliced` and `sparse` segments, which `GetSeq_id` reads, are refused.)
fn query_row_id(align: &Align) -> Result<&ParsedSeqId, ProcessError> {
    let seqalign_exception = |line: u32, kind: &'static str, message: &str| {
        ProcessError::Exception(NcbiException::new(
            ThrowSite {
                file: "src/objects/seqalign/Seq_align.cpp",
                line,
                kind,
                function: None,
            },
            message.to_string(),
        ))
    };
    match &align.segs {
        Segs::Denseg(seg) => seg.ids.first().ok_or_else(|| {
            ProcessError::Unsupported(
                "a Dense-seg without ids (NCBI's CSeq_align::GetSeq_id throws) is not supported by LOSAT's rpsbproc".to_string(),
            )
        }),
        Segs::Std(std) => match &std.query_id {
            QueryRowId::Id(id) => Ok(id),
            QueryRowId::Location(name) => Err(ProcessError::Unsupported(format!(
                "a Std-seg without ids whose first location is a {name} is not supported by LOSAT's rpsbproc"
            ))),
            QueryRowId::Missing => Err(seqalign_exception(
                396,
                "CSeqalignException::eInvalidRowNumber",
                "CSeq_align::GetSeq_id(): can not get seq-id for the row requested.",
            )),
        },
        Segs::Other("packed") => Err(seqalign_exception(
            391,
            "CSeqalignException::eUnsupported",
            "CSeq_align::GetSeq_id() currently does not handle this type of alignment.",
        )),
        Segs::Other(name) => Err(ProcessError::Unsupported(format!(
            "a Seq-align with {name} segments is not supported by LOSAT's rpsbproc"
        ))),
    }
}

/// `ProcessBlastResults`: each query with the alignments that follow in the list while their
/// first id matches the query's local id (else its first id).
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
// each query are not printed. `query_row_id` is `GetSeq_id(0)` of a `Dense-seg`, of the
// `Std-seg`s of a translated search and of a `Packed-seg` (NCBI's exception); the other
// segments are refused there.)
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
        let (accession, title, length, residues) = fill_sequence(query)?;
        let mut dom_query = DomQuery {
            accession,
            title,
            length,
            is_na: query.mol != 3,
            residues,
            aligns: Vec::new(),
            frames: Default::default(),
            translated: None,
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
            let aligned = query_row_id(&aligns[next])?;
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
    use crate::rpsbproc::archive::{StdAlign, StdLoc, StdPair};
    use crate::rpsbproc::data::DataFiles;
    use crate::seqid::model::{Dbtag, TextseqId};

    /// The data files of the fixtures (`tests/fixtures/rpsbproc/data`).
    fn fixture_data() -> CddData {
        let file = |name: &str| {
            format!(
                "{}/tests/fixtures/rpsbproc/data/{name}",
                env!("CARGO_MANIFEST_DIR")
            )
        };
        CddData::load(&DataFiles {
            cddids: file("cddid.tbl"),
            spfeats: file("cddannot.dat"),
            genfeats: file("cddannot_generic.dat"),
            bscores: file("bitscore_specific.txt"),
            cdtrack: file("cdtrack.txt"),
            famlinks: file("family_superfamily_links"),
        })
        .unwrap()
    }

    // The refusals of a nucleotide query's `Std-seg`s (cdalignproc_biodata.cpp:195-232,
    // Seq_align.cpp:339-358 and rpsbproc_main.cpp:1809,1836, quoted at `process_cd_query` and
    // `query_row_id`): a first segment with fewer than 2 ids is read past its end (PS16c oracle
    // `noids1`: segmentation fault in both 0.5.1 builds); a first location of a variant whose
    // id LOSAT does not look up; a kept hit whose range starts or ends at 2^31 - 1, whose
    // DOMAINS line overflows `int` (review oracle `plusmax`), whatever the strand. A hit whose
    // PSSM-ID is 0 is dropped before it is classified or printed.
    #[test]
    fn refusals_of_translated_hits() {
        let data = fixture_data();
        let local = ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()));
        let cdd = |tag: i32| {
            Some(ParsedSeqId::General(Dbtag {
                db: b"CDD".to_vec(),
                tag: ObjectId::Id(tag),
            }))
        };
        let entries = [SeqEntry::Seq(Query {
            ids: vec![local.clone()],
            mol: 1,
            length: Some(6),
            seq_data: Some(("iupacna", b"ACGTAC".to_vec())),
            ..Query::default()
        })];
        let run =
            |query_id: QueryRowId, hit_id: Option<ParsedSeqId>, from: u32, to: u32, strand: i32| {
                let align = Align {
                    segs: Segs::Std(Box::new(StdAlign {
                        query_id,
                        hit_id,
                        segs: vec![StdPair {
                            query: StdLoc::Int {
                                from,
                                to,
                                strand: Some(strand),
                            },
                            hit: StdLoc::Int {
                                from: 0,
                                to: 2,
                                strand: None,
                            },
                        }],
                    })),
                    score: Some(10),
                    e_value: Some(1e-5),
                    bit_score: Some(40.5),
                    num_ident: None,
                    score_of_other_type: None,
                };
                match process_blast_results(&entries, &[align], &data) {
                    Ok(queries) => Ok(queries[0].aligns.len()),
                    Err(ProcessError::Unsupported(message)) => Err(message),
                    Err(other) => panic!("{other:?}"),
                }
            };
        let row = || QueryRowId::Id(local.clone());
        let refused =
            |message: &str| Err(format!("{message} is not supported by LOSAT's rpsbproc"));
        assert_eq!(
            run(row(), None, 0, 2, 1),
            refused("a Seq-align whose first Std-seg has fewer than 2 ids (NCBI reads past them)")
        );
        assert_eq!(
            run(QueryRowId::Location("packed-int"), cdd(5), 0, 2, 1),
            refused("a Std-seg without ids whose first location is a packed-int")
        );
        let end = 0x7fff_ffff;
        let overflow = refused(
            "a hit of a nucleotide query that starts or ends at position 2147483647 (NCBI's int overflows)",
        );
        for strand in [0, 1, 3, 255] {
            assert_eq!(run(row(), cdd(5), 0, end, strand), overflow);
        }
        assert_eq!(run(row(), cdd(5), end, end, 1), overflow);
        assert_eq!(run(row(), cdd(0), 0, end, 1), Ok(0));
    }

    // `Get1LetterSeqData` (biodata_core.cpp:199-253, quoted at `one_letter_residues`): the
    // literals of `ext delta` one after another (NCBIstdaa converted), cut to `inst.length`;
    // the `seq-data` first when both are set; a literal in an encoding LOSAT does not convert is
    // refused as a `seq-data` in it is.
    #[test]
    fn residues_of_delta_literals() {
        let query = |seq_data: Option<(&'static str, Vec<u8>)>,
                     delta: Vec<(&'static str, Vec<u8>)>| Query {
            ids: vec![ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()))],
            mol: 3,
            length: Some(5),
            seq_data,
            delta_literals: delta,
            ..Query::default()
        };
        let residues = |query: Query| fill_sequence(&query).map(|(_, _, _, residues)| residues);
        let delta = vec![
            ("ncbieaa", b"MK".to_vec()),
            ("ncbistdaa", vec![1, 25, 27, 0]),
        ];
        assert_eq!(
            residues(query(None, delta.clone())).ok(),
            Some(b"MKA*J".to_vec())
        );
        assert_eq!(
            residues(query(Some(("iupacaa", b"WW".to_vec())), delta)).ok(),
            Some(b"WW".to_vec())
        );
        let Err(ProcessError::Unsupported(message)) = residues(query(
            None,
            vec![("ncbieaa", b"MK".to_vec()), ("ncbi8aa", vec![1])],
        )) else {
            panic!("a literal of NCBI8aa is refused");
        };
        assert_eq!(
            message,
            "a query whose Seq-data of a delta literal is ncbi8aa is not supported by LOSAT's rpsbproc"
        );
    }

    // `SeqDataTo1LtrSeq` of the nucleotide encodings (biodata_core.cpp:149-197 and the
    // `CSeqportUtil` tables, quoted at `one_letter`): NCBI2na 4 bases per byte from the high
    // bits, NCBI4na 2 (code 0 is `N`), NCBI8na 1 (codes from 16 the byte 0xff), IUPACna as it
    // is; cut to `inst.length`; the literals of a delta of the kind of the first one only
    // (`litType`); NCBIpna throws `CSeqportUtil`'s `runtime_error`.
    #[test]
    fn residues_of_nucleotide_encodings() {
        let query = |length: u32,
                     seq_data: Option<(&'static str, Vec<u8>)>,
                     delta: Vec<(&'static str, Vec<u8>)>| Query {
            ids: vec![ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()))],
            mol: 1,
            length: Some(length),
            seq_data,
            delta_literals: delta,
            ..Query::default()
        };
        let residues = |query: Query| fill_sequence(&query).map(|(_, _, _, residues)| residues);
        // 0x1B = 00 01 10 11, 0xE4 = 11 10 01 00.
        assert_eq!(
            residues(query(7, Some(("ncbi2na", vec![0x1b, 0xe4])), vec![])).ok(),
            Some(b"ACGTTGC".to_vec())
        );
        assert_eq!(
            residues(query(4, Some(("ncbi4na", vec![0x01, 0x2f])), vec![])).ok(),
            Some(b"NACN".to_vec())
        );
        assert_eq!(
            residues(query(3, Some(("ncbi8na", vec![4, 15, 16])), vec![])).ok(),
            Some(vec![b'G', b'N', 0xff])
        );
        assert_eq!(
            residues(query(9, Some(("iupacna", b"acgTUx".to_vec())), vec![])).ok(),
            Some(b"acgTUx".to_vec())
        );
        let delta = vec![
            ("iupacna", b"AC".to_vec()),
            ("ncbieaa", b"MK".to_vec()),
            ("ncbi2na", vec![0xff]),
        ];
        assert_eq!(
            residues(query(10, None, delta)).ok(),
            Some(b"ACTTTT".to_vec())
        );
        let Err(ProcessError::Std(what)) =
            residues(query(4, Some(("ncbipna", vec![0; 5])), vec![]))
        else {
            panic!("NCBIpna throws");
        };
        assert_eq!(what, "Requested conversion not implemented");
    }

    // serial/exception.cpp:86-150 (quoted at `choice_selection`); NCBI's rpsbproc printed the first
    // one for a profile database with `lcl|` ids (PS16a probe 2).
    #[test]
    fn choice_selection_messages() {
        let message = |error| match error {
            ProcessError::Exception(exception) => exception.get_msg().to_string(),
            ProcessError::Unsupported(message) | ProcessError::Std(message) => {
                panic!("{message}")
            }
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
