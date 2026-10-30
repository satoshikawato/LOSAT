//! The definition lines of one database sequence (`Blast-def-line-set`) and the choices that
//! every report makes from them: the Seq-align id (`alignment_id`), the definition line of the
//! subject Bioseq (`bioseq_defline`) and of the shown id (`display_defline`), and the joins of
//! the formats that show every definition line (tabular `sall*`, `salltitles`, `staxids`, the
//! XML `Hit_def`, the lines of an outfmt 0 alignment block).
//!
//! NCBI's formatters read the definition lines through the object manager: the Seq-align id
//! is resolved to the sequence again, whose Bioseq holds the ids and title of one definition
//! line (the first that holds the id) and the whole set as a user object
//! (`CSeqDB::ExtractBlastDefline`). LOSAT knows the sequence and keeps the set; these
//! functions give the parts NCBI takes from each.

use std::borrow::Cow;
use std::collections::BTreeSet;

use super::label::{bare, content, fasta, first_gi, format_label, number, seqid_list_string};
use super::model::{Dbtag, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId, PdbSeqId, TextseqId};
use super::rank::{best_rank_score, find_best, text_score};

/// One `Blast-def-line`, every member decoded (`seqid::ber::decode_full_set`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BlastDefline {
    /// `title`, after NCBI's repair of the stored bytes; `None` when the member is absent
    /// (makeblastdb writes an empty title instead).
    pub title: Option<Vec<u8>>,
    /// The number of stored title bytes outside ' ' ..= '~' that read as '#' (NCBI writes an
    /// error message for each time it reads the definition line).
    pub replaced_title_bytes: usize,
    /// `seqid`, in the stored order.
    pub ids: Vec<ParsedSeqId>,
    /// The number of bytes of the Seq-id strings that read as '#'.
    pub replaced_id_bytes: usize,
    /// `taxid`, as `seqid::ber::taxid_value` reads it.
    pub taxid: Option<i32>,
    pub memberships: Option<Vec<i32>>,
    /// `links`: further taxids of the definition line (`GetTaxIds`).
    pub links: Option<Vec<i32>>,
    pub other_info: Option<Vec<i32>>,
}

/// The `Blast-def-line-set` of one sequence.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DeflineSet {
    pub deflines: Vec<BlastDefline>,
}

/// `CSeq_id::E_SIC`, the result of `CSeq_id::Compare`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sic {
    /// The ids are of types that cannot be compared (`e_DIFF`).
    Diff,
    /// Comparable ids that do not match (`e_NO`).
    No,
    /// The ids match (`e_YES`).
    Yes,
}

/// The definition line whose id an NCBI formatter shows (`display_defline`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DisplayDefline {
    /// The index of the definition line.
    pub index: usize,
    /// Whether it holds the Seq-align id. If none does, NCBI shows the id of the last
    /// definition line and no GI (`gi` is then 0).
    pub found: bool,
    /// The `taxid` of the last definition line up to `index` that has one, 0 for none
    /// (`GetDisplayIds`' `taxid`: it need not be the shown line's).
    pub taxid: i32,
}

/// The kind of id of a tabular field (`CBlastTabularInfo::ESeqIdType`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TabularId {
    /// `sseqid` (`eFullId`).
    Full,
    /// `sacc` (`eAccession`).
    Accession,
    /// `saccver` (`eAccVersion`).
    AccessionVersion,
    /// `sgi` (`eGi`).
    Gi,
}

impl BlastDefline {
    // NCBI reference (598d8ae6): c++/src/objects/blastdb/Blast_def_line.cpp:97-125
    // ```c++
    // CBlast_def_line::GetTaxIds() const
    // {
    //     TTaxIds retval;                 // set<TTaxId>, initially empty
    //
    //     // If there's a 'taxid' value, add it to the result set.
    //     if (CanGetTaxid()) {
    //         TTaxid taxid = GetTaxid();
    //         retval.insert(taxid);
    //     }
    //     // If there are any 'links' values, add them to the result set.
    //     if (IsSetLinks()) {
    //         TLinks taxids = GetLinks();  // see ASN.1 spec comment
    // #ifdef NCBI_STRICT_TAX_ID
    //         ITERATE(TLinks, it, taxids) retval.insert(TAX_ID_FROM(int, *it));
    // #else
    //         retval.insert(taxids.begin(), taxids.end());
    // #endif
    //     }
    //
    //     if(retval.size() > 1) {
    //     	retval.erase(ZERO_TAX_ID);
    //     }
    // ```
    /// `CBlast_def_line::GetTaxIds`: the taxid and the links, without 0 unless it is alone.
    pub fn tax_ids(&self) -> BTreeSet<i32> {
        let mut ids: BTreeSet<i32> = self.taxid.into_iter().collect();
        if let Some(links) = &self.links {
            ids.extend(links.iter().copied());
        }
        if ids.len() > 1 {
            ids.remove(&0);
        }
        ids
    }

    /// The first GI of the definition line, 0 for none (`FindGi`).
    pub fn gi(&self) -> i64 {
        first_gi(&self.ids)
    }

    /// The id of the definition line that NCBI's displays choose: the first `WorstRank`-best.
    pub fn display_id(&self) -> Option<&ParsedSeqId> {
        find_best(&self.ids, text_score)
    }
}

impl DeflineSet {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1844-1864
    // ```c++
    // list< CRef<CSeq_id> > CSeqDBVol::GetSeqIDs(int  oid) const
    // {
    //     list< CRef< CSeq_id > > seqids;
    //
    //     CRef<CBlast_def_line_set> defline_set =
    //         x_GetFilteredHeader(oid, NULL);
    //
    //     if ((! defline_set.Empty()) && defline_set->CanGet()) {
    //         ITERATE(list< CRef<CBlast_def_line> >, defline, defline_set->Get()) {
    //             if (! (*defline)->CanGetSeqid()) {
    //                 continue;
    //             }
    //
    //             ITERATE(list< CRef<CSeq_id> >, seqid, (*defline)->GetSeqid()) {
    //                 seqids.push_back(*seqid);
    //             }
    //         }
    //     }
    //
    //     return seqids;
    // }
    // ```
    // `x_GetFilteredHeader` drops definition lines only for a membership bit or a GI list,
    // which LOSAT's searches do not use.
    /// Every Seq-id of every definition line, in order (`CSeqDB::GetSeqIDs`).
    pub fn ids(&self) -> impl Iterator<Item = &ParsedSeqId> {
        self.deflines.iter().flat_map(|defline| defline.ids.iter())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqinfosrc_aux.cpp:51-69
    // ```c++
    // void GetSequenceLengthAndId(const blast::IBlastSeqInfoSrc * seqinfo_src,
    //                             int                      oid,
    //                             CRef<CSeq_id>          & seqid,
    //                             TSeqPos                * length)
    // {
    //     _ASSERT(length);
    //     list<CRef<CSeq_id> > seqid_list = seqinfo_src->GetId(oid);
    //
    //     CRef<CSeq_id> id = FindBestChoice(seqid_list, CSeq_id::BestRank);
    //     if (id.NotEmpty()) {
    //         seqid.Reset(new CSeq_id);
    // 	if (id->IsGi())
    //         	SerialAssign(*seqid, *id);
    // 	else
    //         	SerialAssign(*seqid, *(seqid_list.front()));
    //     }
    // ```
    /// The subject Seq-id of the sequence's Seq-aligns (outfmt 8, 9, 11, 12 print it; NCBI's
    /// formatters find the subject by it): the first GI of any definition line, else the first
    /// id of the first definition line. `None` for a set without ids.
    pub fn alignment_id(&self) -> Option<&ParsedSeqId> {
        let best = find_best(self.ids(), best_rank_score)?;
        if best.gi().is_some() {
            Some(best)
        } else {
            self.ids().next()
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1328-1359
    // ```c++
    //     if ((target_gi != ZERO_GI) || target_seq_id) {
    //         defline_set.Reset(new CBlast_def_line_set);
    //
    //         CRef<const CSeq_id > seqid;
    //         CRef<const CSeq_id > seqid_lcl;
    //         if (target_gi != ZERO_GI) {
    //             seqid.Reset(new CSeq_id(CSeq_id::e_Gi, target_gi));
    //             seqid_lcl.Reset(new CSeq_id(CSeq_id::e_Local, NStr::NumericToString(target_gi)));
    //         }
    //         else {
    //             seqid.Reset(target_seq_id);
    //             seqid_lcl.Reset(new CSeq_id(CSeq_id::e_Local, seqid->GetSeqIdString(true)));
    //         }
    //
    //         CRef<CBlast_def_line> filt_dl;
    //
    //         ITERATE(TDeflines, iter, orig_deflines->Get()) {
    //             if ((seqid.NotEmpty() && s_SeqDB_SeqIdIn((**iter).GetSeqid(), *seqid)) ||
    //                 (seqid_lcl.NotEmpty() && s_SeqDB_SeqIdIn((**iter).GetSeqid(), *seqid_lcl))) {
    //                 filt_dl = *iter;
    //                 break;
    //             }
    //         }
    //
    //         if (filt_dl.Empty()) {
    //             NCBI_THROW(CSeqDBException, eArgErr,
    //                        "Error: oid headers do not contain target gi/seq_id.");
    // ```
    // The data loader asks for the Bioseq of the Seq-align id (cached_sequence.cpp:58):
    // `target` is the `alignment_id` of the same set.
    /// The definition line of the subject Bioseq that NCBI's formatters read for `target` (the
    /// Seq-align id): the first that holds it (`s_SeqDB_SeqIdIn`, or a local id of its content
    /// label). Its ids and title are the Bioseq's (tabular `sseqid`, `sacc`, ...). `None` when
    /// no definition line holds it (NCBI throws).
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2898-2916
    // ```c++
    // CSeq_id& CSeq_id::Set(E_Choice           the_type,
    //                       const CTempString& acc_in,
    //                       ...
    // {
    //     CTempString  acc       = NStr::TruncateSpaces_Unsafe(acc_in,
    //                                                          NStr::eTrunc_Both);
    //     ...
    //     case e_Local:
    //         SetLocal().SetStrOrId(acc);
    // ```
    // `TruncateSpaces` drops the bytes that C's `isspace` takes; `SetStrOrId` gives an
    // integer tag for a label that reads as one, which `compare` treats as the same tag.
    pub fn bioseq_defline(&self, target: &ParsedSeqId) -> Option<usize> {
        let mut label = Vec::new();
        content(target, true, &mut label);
        let is_space = |c: &u8| matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r');
        let start = label
            .iter()
            .position(|c| !is_space(c))
            .unwrap_or(label.len());
        let end = label
            .iter()
            .rposition(|c| !is_space(c))
            .map_or(start, |end| end + 1);
        let local = ParsedSeqId::Local(ObjectId::Str(label[start..end].to_vec()));
        self.deflines
            .iter()
            .position(|defline| seqid_in(&defline.ids, target) || seqid_in(&defline.ids, &local))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4741-4785
    // ```c++
    //     if(gi) *gi = ZERO_GI;
    //     if(taxid) *taxid = ZERO_TAX_ID;
    //     ...
    //         for(list< CRef< CBlast_def_line > >::const_iterator iter = bdl.begin();
    //             iter != bdl.end(); iter++){
    //             const CBioseq::TId* cur_id = &((*iter)->GetSeqid());
    //             TGi cur_gi =  FindGi(*cur_id);
    //             wid = FindBestChoice(*cur_id, CSeq_id::WorstRank);
    //             string curSeqID = GetLabel(wid,true);//uses GetSeqIdString(true)
    //             if (taxid && (*iter)->IsSetTaxid() && (*iter)->CanGetTaxid()){
    //                 *taxid = (*iter)->GetTaxid();
    //             }
    //             ...
    //             } else {
    //                 ITERATE(CBioseq::TId, iter_id, *cur_id) {
    //                     if ((*iter_id)->Match(aln_id)
    //                       || (aln_id.IsGeneral() && aln_id.GetGeneral().CanGetDb() &&
    //                          (*iter_id)->IsGeneral() && (*iter_id)->GetGeneral().CanGetDb() &&
    //                          aln_id.GetGeneral().GetDb() == (*iter_id)->GetGeneral().GetDb())) {
    //                         found = true;
    //                     }
    //                 }
    //             }
    //             if(found){
    //                 if(gi) *gi = cur_gi;
    //                 if(textSeqID) *textSeqID = curSeqID;
    //                 break;
    //             }
    //         }
    //     }
    //
    //     return wid;
    // ```
    /// The definition line whose `WorstRank`-best id NCBI shows for the Seq-align id `aln_id`
    /// (`CAlignFormatUtil::GetDisplayIds`: outfmt 0 descriptions): the first that holds
    /// `aln_id` or a `general` id of the same database, else the last one. `None` for a set
    /// without definition lines.
    pub fn display_defline(&self, aln_id: &ParsedSeqId) -> Option<DisplayDefline> {
        let last = self.deflines.len().checked_sub(1)?;
        let found = self.deflines.iter().position(|defline| {
            defline.ids.iter().any(|id| {
                compare(id, aln_id) == Sic::Yes
                    || matches!((id, aln_id),
                        (ParsedSeqId::General(a), ParsedSeqId::General(b)) if a.db == b.db)
            })
        });
        let index = found.unwrap_or(last);
        let taxid = self.deflines[..=index]
            .iter()
            .rev()
            .find_map(|defline| defline.taxid)
            .unwrap_or(0);
        Some(DisplayDefline {
            index,
            found: found.is_some(),
            taxid,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:590-603
    // ```c++
    //     if (bdlRef.NotEmpty() && bdlRef->CanGet() && bdlRef->IsSet() && !bdlRef->Get().empty()){
    //     	ITERATE(CBlast_def_line_set::Tdata, itr, bdlRef->Get()) {
    //             if((*itr)->IsSetTaxid()) {
    //                     if((*itr)->GetTaxid() != ZERO_TAX_ID) {
    //                     	m_SubjectTaxId = (*itr)->GetTaxid();
    //                     	break;
    //                     }
    //             }
    //     	}
    //     }
    //
    //     if(m_SubjectTaxId == ZERO_TAX_ID) {
    //           m_SubjectTaxId = sequence::GetTaxId(handle);
    //     }
    // ```
    // The fall-back finds nothing: the Bioseq's taxonomy descriptors are those of its
    // definition line's taxid when it is above 0 (`x_GetTaxonomy`, seqdbvol.cpp:1131-1137),
    // and then the loop has found a non-zero taxid already.
    /// The `staxid` of the set: the first non-zero `taxid` (`x_SetTaxInfo`), 0 for none.
    pub fn first_taxid(&self) -> i32 {
        self.deflines
            .iter()
            .filter_map(|defline| defline.taxid)
            .find(|&taxid| taxid != 0)
            .unwrap_or(0)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:647-670
    // ```c++
    //     if (bdlRef.NotEmpty() && bdlRef->CanGet() && bdlRef->IsSet() && !bdlRef->Get().empty()){
    //     	ITERATE(CBlast_def_line_set::Tdata, itr, bdlRef->Get()) {
    //     		CBlast_def_line::TTaxIds t = (*itr)->GetTaxIds();
    //     		m_SubjectTaxIds.insert(t.begin(), t.end());
    //     	}
    //     }
    //
    //     if(m_SubjectTaxIds.empty()) {
    //             CSeqdesc_CI desc_s(handle, CSeqdesc::e_Source);
    //             for (;desc_s; ++desc_s) {
    //                      TTaxId t = desc_s->GetSource().GetOrg().GetTaxId();
    //                      if(t != ZERO_TAX_ID) {
    //                     	 m_SubjectTaxIds.insert(t);
    //                      }
    //             }
    //
    //             CSeqdesc_CI desc(handle, CSeqdesc::e_Org);
    //             for (; desc; ++desc) {
    //                     TTaxId t= desc->GetOrg().GetTaxId();
    //                     if(t != ZERO_TAX_ID) {
    //                     	m_SubjectTaxIds.insert(t);
    //                     }
    //             }
    //     }
    // ```
    // The fall-back finds nothing for the reason given at `first_taxid`: a definition line
    // with a taxid above 0 makes `GetTaxIds` non-empty. `m_SubjectTaxIds` is a `set`.
    /// The `staxids` of the set: the union of every definition line's `GetTaxIds`
    /// (`x_SetTaxInfoAll`), in ascending order.
    pub fn tax_ids(&self) -> BTreeSet<i32> {
        let mut ids = BTreeSet::new();
        for defline in &self.deflines {
            ids.extend(defline.tax_ids());
        }
        ids
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:544-555
    // ```c++
    //     if (bdlRef.NotEmpty() && bdlRef->CanGet() && bdlRef->IsSet() && !bdlRef->Get().empty()){
    //         vector< CConstRef<CSeq_id> > original_seqids;
    //
    //         ITERATE(CBlast_def_line_set::Tdata, itr, bdlRef->Get()) {
    //             original_seqids.clear();
    //             ITERATE(CBlast_def_line::TSeqid, id, (*itr)->GetSeqid()) {
    //                 original_seqids.push_back(*id);
    //             }
    //             list<CRef<objects::CSeq_id> > next_seqid_list;
    //             // Next call replaces BL_ORD_ID if found.
    //             CShowBlastDefline::GetSeqIdList(bh,original_seqids,next_seqid_list);
    //             m_SubjectIds.push_back(next_seqid_list);
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:224-231
    // ```c++
    // void CBlastTabularInfo::x_PrintSubjectAllSeqIds(void)
    // {
    //     ITERATE(vector<list<CRef<CSeq_id> > >, iter, m_SubjectIds) {
    //         if (iter != m_SubjectIds.begin())
    //             m_Ostream << ";";
    //         m_Ostream << s_GetSeqIdListString(*iter, eFullId);
    //     }
    // }
    // ```
    // `x_PrintSubjectAllGis` and `x_PrintSubjectAllAccessions` (tabular.cpp:238-264) are the
    // same loop with `eGi` and `eAccession`.
    /// The tabular field of `kind` for each definition line, `;`-joined (`sallseqid`,
    /// `sallacc`, `sallgi`): `tabular_id` of the line's ids after `seq_id_list` with
    /// `title_word` (the subject Bioseq's, the same for every line).
    pub fn tabular_all(&self, kind: TabularId, title_word: &[u8], out: &mut Vec<u8>) {
        for (i, defline) in self.deflines.iter().enumerate() {
            if i > 0 {
                out.push(b';');
            }
            tabular_id(&seq_id_list(&defline.ids, title_word), kind, out);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:387-440
    // ```c++
    // void CBlastTabularInfo::x_PrintSubjectAllTitles ()
    // {
    // 	if(m_SubjectDefline.NotEmpty() && m_SubjectDefline->CanGet() &&
    // 	   m_SubjectDefline->IsSet() && !m_SubjectDefline->Get().empty())
    // 	{
    // 		const list<CRef<CBlast_def_line> > & defline = m_SubjectDefline->Get();
    // 		list<CRef<CBlast_def_line> >::const_iterator iter = defline.begin();
    // 		for(; iter != defline.end(); ++iter)
    // 		{
    // 			if (iter != defline.begin())
    // 						 m_Ostream << "<>";
    //
    // 			if((*iter)->IsSetTitle())
    // 			{
    // 				if((*iter)->GetTitle().empty())
    // 					m_Ostream << NA;
    // 				else
    // 					m_Ostream << (*iter)->GetTitle();
    // 			}
    // 			else
    // 				m_Ostream << NA;
    // 		}
    // 	}
    // 	else
    // 		m_Ostream << NA;
    //
    // }
    //
    // void CBlastTabularInfo::x_PrintSubjectTitle()
    // ...
    // 			if(defline.front()->IsSetTitle())
    // 			{
    // 				if(defline.front()->GetTitle().empty())
    // 					m_Ostream << NA;
    // 				else
    // 					m_Ostream << defline.front()->GetTitle();
    // 			}
    // 			else
    // 				m_Ostream << NA;
    // ```
    /// `stitle`: the first title, `N/A` for an empty or absent one or no definition line.
    pub fn first_title(&self, out: &mut Vec<u8>) {
        tabular_title(self.deflines.first(), out);
    }

    /// `salltitles`: every title as `stitle` prints it, `<>`-joined.
    pub fn all_titles(&self, out: &mut Vec<u8>) {
        if self.deflines.is_empty() {
            out.extend_from_slice(NOT_AVAILABLE);
        }
        for (i, defline) in self.deflines.iter().enumerate() {
            if i > 0 {
                out.extend_from_slice(b"<>");
            }
            tabular_title(Some(defline), out);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:266-304
    // ```c++
    //         bdlRef->PutTargetGiFirst(this_gi_first);
    //         const list< CRef< CBlast_def_line > >& bdl = bdlRef->Get();
    //         bool is_first = true;
    //         ITERATE(list<CRef<CBlast_def_line> >, iter, bdl) {
    //             const CBioseq::TId& cur_id = (*iter)->GetSeqid();
    //             TGi cur_gi =  FindGi(cur_id);
    //             TGi gi_in_use_this_gi = ZERO_GI;
    //             ...
    //             if(use_this_gi.empty() || gi_in_use_this_gi > ZERO_GI) {
    //                 if (is_first)
    //                     seqid = GetSeqIdListString(cur_id, show_gi);
    //
    //                 if((*iter)->IsSetTitle()){
    //                     if(is_first){
    //                         defline = (*iter)->GetTitle();
    //                     } else {
    //                         string concat_acc;
    //                         CConstRef<CSeq_id> wid =
    //                             FindBestChoice(cur_id, CSeq_id::WorstRank);
    //                         wid->GetLabel(&concat_acc, CSeq_id::eFasta, 0);
    //                         if( show_gi && cur_gi > ZERO_GI){
    //                             defline = defline + " >" + "gi|" +
    //                                 NStr::NumericToString(cur_gi) + "|" +
    //                                 concat_acc + " " + (*iter)->GetTitle();
    //                         } else {
    //                             defline = defline + " >" + concat_acc + " " +
    //                                 (*iter)->GetTitle();
    //                         }
    //                     }
    //                     is_first = false;
    //                 }
    //             }
    //         }
    // ```
    // `PutTargetGiFirst(-1)` (the default `this_gi_first`) leaves the set as it is
    // (Blast_def_line_set.cpp:212-237).
    /// `CShowBlastDefline::GetBioseqHandleDeflineAndId` (outfmt 5 `Hit_id` and `Hit_def`):
    /// the id list (`seqid_list_string`) of the first definition line with a title (of the
    /// last one when none has a title), and that title followed by ` >[gi|N|]<fasta> <title>`
    /// for each later definition line with a title.
    pub fn hit_id_and_def(&self, show_gi: bool) -> (Vec<u8>, Vec<u8>) {
        let (mut id, mut title) = (Vec::new(), Vec::new());
        let mut is_first = true;
        for defline in &self.deflines {
            if is_first {
                id.clear();
                seqid_list_string(&defline.ids, show_gi, &mut id);
            }
            let Some(text) = &defline.title else {
                continue;
            };
            if !is_first {
                title.extend_from_slice(b" >");
                let gi = defline.gi();
                if show_gi && gi > 0 {
                    title.extend_from_slice(b"gi|");
                    number(&mut title, gi);
                    title.push(b'|');
                }
                if let Some(best) = defline.display_id() {
                    fasta(best, &mut title);
                }
                title.push(b' ');
            }
            title.extend_from_slice(text);
            is_first = false;
        }
        (id, title)
    }
}

const NOT_AVAILABLE: &[u8] = b"N/A";

fn tabular_title(defline: Option<&BlastDefline>, out: &mut Vec<u8>) {
    match defline.and_then(|defline| defline.title.as_deref()) {
        Some(title) if !title.is_empty() => out.extend_from_slice(title),
        _ => out.extend_from_slice(NOT_AVAILABLE),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:165-197
// ```c++
// static string
// s_GetSeqIdListString(const list<CRef<CSeq_id> >& id,
//                      CBlastTabularInfo::ESeqIdType id_type)
// {
//     string id_str = NcbiEmptyString;
//
//     switch (id_type) {
//     case CBlastTabularInfo::eFullId:
//         id_str = CShowBlastDefline::GetSeqIdListString(id, true);
//         break;
//     case CBlastTabularInfo::eAccession:
//     {
//         CConstRef<CSeq_id> accid = FindBestChoice(id, CSeq_id::Score);
//         accid->GetLabel(&id_str, CSeq_id::eContent, 0);
//         break;
//     }
//     case CBlastTabularInfo::eAccVersion:
//     {
//         CConstRef<CSeq_id> accid = FindBestChoice(id, CSeq_id::Score);
//         accid->GetLabel(&id_str, CSeq_id::eContent, CSeq_id::fLabel_Version);
//         break;
//     }
//     case CBlastTabularInfo::eGi:
//         id_str = NStr::NumericToString(FindGi(id));
//         break;
//     default: break;
//     }
//
//     if (id_str == NcbiEmptyString)
//         id_str = "Unknown";
//
//     return id_str;
// }
// ```
// NCBI dereferences the best id of an empty list (accession kinds); every definition line
// that makeblastdb writes has an id, and LOSAT writes `Unknown` for none.
/// The tabular id of `kind` of the id list `ids` (`s_GetSeqIdListString`): `Unknown` for an
/// empty result. The subject fields (`sseqid`, `sacc`, `saccver`, `sgi`) take the ids of the
/// subject Bioseq's definition line (`bioseq_defline`) after `seq_id_list`.
pub fn tabular_id(ids: &[ParsedSeqId], kind: TabularId, out: &mut Vec<u8>) {
    let start = out.len();
    match kind {
        TabularId::Full => seqid_list_string(ids, true, out),
        TabularId::Accession | TabularId::AccessionVersion => {
            if let Some(best) = find_best(ids, text_score) {
                content(best, kind == TabularId::AccessionVersion, out);
            }
        }
        TabularId::Gi => number(out, first_gi(ids)),
    }
    if out.len() == start {
        out.extend_from_slice(b"Unknown");
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:212-243
// ```c++
// void
// CShowBlastDefline::GetSeqIdList(const objects::CBioseq_Handle& bh,
//                                 vector< CConstRef<CSeq_id> > &original_seqids,
//                                 list<CRef<objects::CSeq_id> >& ids)
// {
//     ids.clear();
//     ITERATE(vector< CConstRef<CSeq_id> >, itr, original_seqids) {
//         CRef<CSeq_id> next_seqid(new CSeq_id());
//         string id_token = NcbiEmptyString;
//
//         if (((*itr)->IsGeneral() &&
//             (*itr)->AsFastaString().find("gnl|BL_ORD_ID")
//             != string::npos) ||
// 		(*itr)->AsFastaString().find("lcl|Subject_") != string::npos) {
//             vector<string> title_tokens;
//             string defline = sequence::CDeflineGenerator().GenerateDefline(bh);
//             if (defline != NcbiEmptyString) {
//                 id_token =
//                     NStr::Split(defline, " ", title_tokens)[0];
//             }
//         }
//         if (id_token != NcbiEmptyString) {
//             // Create a new local id with a label containing the extracted
//             // token and save it in the next_seqid instead of the original
//             // id.
//             CObject_id* obj_id = new CObject_id();
//             obj_id->SetStr(id_token);
//             next_seqid->SetLocal(*obj_id);
//         } else {
//             next_seqid->Assign(**itr);
//         }
//         ids.push_back(next_seqid);
//     }
// }
// ```
// The tabular formats call it for the subject Bioseq's ids (`SetSubjectId`, tabular.cpp:
// 522-532; its `s_ReplaceLocalId` keeps a local id's own text for a `-db` search, where
// `m_ParseSubjectDefline` is false) and for every definition line (`x_SetSubjectIds`).
/// `CShowBlastDefline::GetSeqIdList`: `ids` with each `general` id whose FASTA form holds
/// `gnl|BL_ORD_ID`, and each id whose FASTA form holds `lcl|Subject_`, replaced by
/// `local str "<title_word>"`; with an empty `title_word` the ids stay. `title_word` is the
/// first word (the bytes before the first space) of the subject Bioseq's generated title,
/// `CDeflineGenerator().GenerateDefline(bh)`, which the caller makes: the Bioseq's title is
/// its definition line's (seqdbvol.cpp:1445-1463; `unnamed protein product` for an empty
/// protein title, empty for an empty nucleotide one), and its Seq-ids can add a prefix
/// (`TPA: ` for a `tpg`, `tpe` or `tpd` id) or a generated title (a `patent` or `pdb` id
/// with an empty title).
pub fn seq_id_list<'a>(ids: &'a [ParsedSeqId], title_word: &[u8]) -> Cow<'a, [ParsedSeqId]> {
    let replaced = |id: &ParsedSeqId| {
        let mut fasta_form = Vec::new();
        fasta(id, &mut fasta_form);
        (matches!(id, ParsedSeqId::General(_)) && contains(&fasta_form, b"gnl|BL_ORD_ID"))
            || contains(&fasta_form, b"lcl|Subject_")
    };
    if title_word.is_empty() || !ids.iter().any(replaced) {
        return Cow::Borrowed(ids);
    }
    Cow::Owned(
        ids.iter()
            .map(|id| {
                if replaced(id) {
                    ParsedSeqId::Local(ObjectId::Str(title_word.to_vec()))
                } else {
                    id.clone()
                }
            })
            .collect(),
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2417-2437
// ```c++
//                     if(m_AlignOption&eShowGi && alnDispParams->gi > ZERO_GI &&
//                        !alnDispParams->seqID->IsGi()){
//                         out<<"gi|"<<alnDispParams->gi<<"|";
//                     }
//                     if(!(alnDispParams->seqID->AsFastaString().find("gnl|BL_ORD_ID") != string::npos) ||
// 			alnDispParams->seqID->AsFastaString().find("lcl|Subject_") != string::npos){
//                     	if (strncmp(alnDispParams->seqID->AsFastaString().c_str(), "lcl|", 4) == 0) {
//                             out << alnDispParams->label;
//                     	}
//                     	else {
//                             if (m_UseLongSeqIds ||
//                                 ((m_AlignOption & eShowGi) &&
//                                  alnDispParams->gi > ZERO_GI)) {
//                                 alnDispParams->seqID->WriteAsFasta(out);
//                             }
//                             else {
//                                 out << CAlignFormatUtil::GetBareId(
//                                                        *alnDispParams->seqID);
//                             }
//                     	}
//                     }
// ```
// `seqID` is the `WorstRank`-best id of the definition line, `gi` its first GI, `label`
// `CAlignFormatUtil::GetLabel(seqID, CSeq_id::eContent)` (the enum is the `bool
// with_version` parameter: true) (showalign.cpp:2205-2209). `m_UseLongSeqIds` is the
// `[BLAST] LONG_SEQID` registry value, false by default.
/// The id that an outfmt 0 alignment block prints for `defline` (after `>` or ` `, before
/// ` <title>`; `CDisplaySeqalign::x_PrintDefLine`): `gi|N|` and the FASTA form with
/// `-show_gis` and a GI, else the bare id (a local id as its label).
pub fn alignment_line_id(defline: &BlastDefline, show_gis: bool, out: &mut Vec<u8>) {
    let Some(id) = defline.display_id() else {
        return;
    };
    let gi = defline.gi();
    if show_gis && gi > 0 && id.gi().is_none() {
        out.extend_from_slice(b"gi|");
        number(out, gi);
        out.push(b'|');
    }
    let mut fasta_form = Vec::new();
    fasta(id, &mut fasta_form);
    if !contains(&fasta_form, b"gnl|BL_ORD_ID") || contains(&fasta_form, b"lcl|Subject_") {
        if fasta_form.starts_with(b"lcl|") {
            format_label(id, true, out);
        } else if show_gis && gi > 0 {
            out.extend_from_slice(&fasta_form);
        } else {
            bare(id, out);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:886-906
// ```c++
//         if(m_Option & eShowGi){
//             if(sdl->gi > ZERO_GI){
//                 line_component = "gi|" + NStr::NumericToString(sdl->gi) + "|";
//                 out << line_component;
//                 line_length += line_component.size();
//             }
//         }
//         if(!sdl->id.Empty()){
//             if(!(sdl->id->AsFastaString().find("gnl|BL_ORD_ID") != string::npos ||
// 		sdl->id->AsFastaString().find("lcl|Subject_") != string::npos)){
//                 string idStr;
//                 if (use_long_seqids || ((m_Option & eShowGi) && !sdl->id->IsGi())) {
//                     idStr = sdl->id->AsFastaString();
//                 }
//                 else {
//                     idStr = CAlignFormatUtil::GetBareId(*sdl->id) + " ";
//                 }
//             	if (strncmp(idStr.c_str(), "lcl|", 4) == 0) {
//             		idStr = sdl->id->AsFastaString().substr(4);
//             	}
//                 out << idStr;
//                 line_length += idStr.size();
//             }
//         }
//         ...
//         line_component = (line_component.empty() ? "" : "  ") + sdl->defline;
// ```
// `sdl->id` and `sdl->gi` are the id and GI of `GetDisplayIds` (`display_defline`; the GI
// is 0 unless the definition line holds the Seq-align id).
/// The start of an outfmt 0 description line (`CShowBlastDefline`, text mode): `gi|N|` with
/// `-show_gis` and a GI, then the id (`GetBareId` and a space; the FASTA form with `-show_gis`
/// unless the id is the GI; a local id without `lcl|`). Returns whether `gi|N|` was written:
/// NCBI then puts two spaces before the title.
pub fn description_id(id: &ParsedSeqId, gi: i64, show_gis: bool, out: &mut Vec<u8>) -> bool {
    let wrote_gi = show_gis && gi > 0;
    if wrote_gi {
        out.extend_from_slice(b"gi|");
        number(out, gi);
        out.push(b'|');
    }
    let mut fasta_form = Vec::new();
    fasta(id, &mut fasta_form);
    if !(contains(&fasta_form, b"gnl|BL_ORD_ID") || contains(&fasta_form, b"lcl|Subject_")) {
        let mut id_str = Vec::new();
        if show_gis && id.gi().is_none() {
            id_str.extend_from_slice(&fasta_form);
        } else {
            bare(id, &mut id_str);
            id_str.push(b' ');
        }
        if id_str.starts_with(b"lcl|") {
            out.extend_from_slice(&fasta_form[4..]);
        } else {
            out.extend_from_slice(&id_str);
        }
    }
    wrote_gi
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    haystack
        .windows(needle.len())
        .any(|window| window == needle)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1031-1052
// ```c++
// s_SeqDB_SeqIdIn(const list< CRef< CSeq_id > > & seqids, const CSeq_id & target)
// {
//     typedef list< CRef<CSeq_id> > TSeqidList;
//
//     ITERATE(TSeqidList, iter, seqids) {
//         CSeq_id::E_SIC rv = (**iter).Compare(target);
//
//         switch(rv) {
//         case CSeq_id::e_YES:
//             return true;
//
//         case CSeq_id::e_NO:
//             return false;
//
//         default:
//             break;
//         }
//     }
//
//     return false;
// }
// ```
/// `s_SeqDB_SeqIdIn`: whether the first id of `ids` that is comparable with `target` matches.
fn seqid_in(ids: &[ParsedSeqId], target: &ParsedSeqId) -> bool {
    for id in ids {
        match compare(id, target) {
            Sic::Yes => return true,
            Sic::No => return false,
            Sic::Diff => {}
        }
    }
    false
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:392-482
// ```c++
// inline bool CanCmpAcc(CSeq_id::E_Choice choice)
// {
//     switch ( choice ) {
//     case CSeq_id::e_Genbank:
//     case CSeq_id::e_Embl:
//     case CSeq_id::e_Ddbj:
//     case CSeq_id::e_Tpg:
//     case CSeq_id::e_Tpe:
//     case CSeq_id::e_Tpd:
//     case CSeq_id::e_Gpipe:
//     case CSeq_id::e_Named_annot_track:
//         return true;
//     default:
//         return false;
//     }
// }
// ...
// CSeq_id::E_SIC CSeq_id::Compare(const CSeq_id& sid2) const
// {
//     if ( Which() != sid2.Which() ) { // Only one case where this will work
//         if (!CanCmpAcc(Which()) || !CanCmpAcc(sid2.Which())) {
//             return e_DIFF;
//         }
//         const CTextseq_id *tsip1 = GetTextseq_Id();
//         ...
//         if ( tsip1->Match(*tsip2) ) // id Textseq_id match
//             return e_YES;
//         else
//             return e_NO;
//     }
//
//     switch ( Which() ) { // Now we only need to know one
//     case e_Local:
//         return GetLocal().Match(sid2.GetLocal()) ? e_YES : e_NO;
//     case e_Gibbsq:
//         return GetGibbsq() == sid2.GetGibbsq() ? e_YES : e_NO;
//     case e_Gibbmt:
//         return GetGibbmt() == sid2.GetGibbmt() ? e_YES : e_NO;
//     case e_Giim:
//         return GetGiim().GetId() == sid2.GetGiim().GetId() ? e_YES : e_NO;
//     ...
//     case e_Patent:
//         return GetPatent().Match(sid2.GetPatent()) ? e_YES : e_NO;
//         (CPatent_seq_id::Match, CId_pat::Match: quoted at `patent_match`)
//     ...
//     case e_General:
//         if ( GetGeneral().Match(sid2.GetGeneral()) ) {
//             return e_YES;
//         }
//         else if ( NStr::CompareNocase(GetGeneral().GetDb(),
//             sid2.GetGeneral().GetDb()) ) {
//                 return e_DIFF;
//         }
//         return e_NO;
//     case e_Gi:
//         return GetGi() == sid2.GetGi() ? e_YES : e_NO;
//     ...
//     case e_Pdb:
//         return GetPdb().Match(sid2.GetPdb()) ? e_YES : e_NO;
//     (every Textseq-id alternative: its Textseq-id's Match)
// ```
/// `CSeq_id::Compare`.
pub fn compare(a: &ParsedSeqId, b: &ParsedSeqId) -> Sic {
    let yes_no = |matches: bool| if matches { Sic::Yes } else { Sic::No };
    if a.alternative() != b.alternative() {
        let comparable = |id: &ParsedSeqId| {
            matches!(
                id,
                ParsedSeqId::Genbank(_)
                    | ParsedSeqId::Embl(_)
                    | ParsedSeqId::Ddbj(_)
                    | ParsedSeqId::Tpg(_)
                    | ParsedSeqId::Tpe(_)
                    | ParsedSeqId::Tpd(_)
                    | ParsedSeqId::Gpipe(_)
                    | ParsedSeqId::NamedAnnotTrack(_)
            )
        };
        if !comparable(a) || !comparable(b) {
            return Sic::Diff;
        }
        return yes_no(textseq_match(a.textseq().unwrap(), b.textseq().unwrap()));
    }
    match (a, b) {
        (ParsedSeqId::Local(x), ParsedSeqId::Local(y)) => yes_no(object_id_match(x, y)),
        (ParsedSeqId::Gibbsq(x), ParsedSeqId::Gibbsq(y))
        | (ParsedSeqId::Gibbmt(x), ParsedSeqId::Gibbmt(y)) => yes_no(x == y),
        (ParsedSeqId::Giim(x), ParsedSeqId::Giim(y)) => yes_no(x.id == y.id),
        (ParsedSeqId::Patent(x), ParsedSeqId::Patent(y)) => yes_no(patent_match(x, y)),
        (ParsedSeqId::General(x), ParsedSeqId::General(y)) => {
            if dbtag_match(x, y) {
                Sic::Yes
            } else if !x.db.eq_ignore_ascii_case(&y.db) {
                Sic::Diff
            } else {
                Sic::No
            }
        }
        (ParsedSeqId::Gi(x), ParsedSeqId::Gi(y)) => yes_no(x == y),
        (ParsedSeqId::Pdb(x), ParsedSeqId::Pdb(y)) => yes_no(pdb_match(x, y)),
        _ => yes_no(textseq_match(a.textseq().unwrap(), b.textseq().unwrap())),
    }
}

// No NCBI counterpart: a necessary condition of `compare`'s `Yes`, for LOSAT's check that no two
// sequences of a database can be taken for one where a report joins following hits
// (`report::db_subject::join_keys`).
/// Adds to `keys` the keys of `id`: two ids that `compare` finds the same (`Sic::Yes`) share at
/// least one key, so ids that share none never do (ids that share one need not). A text Seq-id
/// has its accession and its name in upper case, whatever its alternative (a superset of the
/// alternatives `compare` compares across); a `pdb` id its molecule; the other alternatives their
/// alternative and their value, an object id as its number or its string in upper case
/// (`object_id_match`). The keys are 64-bit hashes of those values.
pub fn match_keys(id: &ParsedSeqId, keys: &mut Vec<u64>) {
    use std::hash::{Hash, Hasher};
    let mut key = |parts: &[&[u8]]| {
        let mut hasher = std::hash::DefaultHasher::new();
        parts.hash(&mut hasher);
        keys.push(hasher.finish());
    };
    let tag_key = |tag: &ObjectId| -> Vec<u8> {
        match object_id_number(tag) {
            Some(number) => format!("n{number}").into_bytes(),
            None => match tag {
                ObjectId::Str(text) => [b"s".as_slice(), &text.to_ascii_uppercase()].concat(),
                ObjectId::Id(_) => unreachable!("an integer tag is a number"),
            },
        }
    };
    match id {
        ParsedSeqId::Local(tag) => key(&[b"local", &tag_key(tag)]),
        ParsedSeqId::Gibbsq(number) => key(&[b"gibbsq", &number.to_le_bytes()]),
        ParsedSeqId::Gibbmt(number) => key(&[b"gibbmt", &number.to_le_bytes()]),
        ParsedSeqId::Giim(giim) => key(&[b"giim", &giim.id.to_le_bytes()]),
        ParsedSeqId::Patent(patent) => {
            let number = match &patent.cit.id {
                PatentNumber::Number(number) | PatentNumber::AppNumber(number) => number,
            };
            key(&[
                b"patent",
                &patent.seqid.to_le_bytes(),
                &patent.cit.country.to_ascii_uppercase(),
                &number.to_ascii_uppercase(),
            ])
        }
        ParsedSeqId::General(dbtag) => key(&[
            b"general",
            &dbtag.db.to_ascii_uppercase(),
            &tag_key(&dbtag.tag),
        ]),
        ParsedSeqId::Gi(gi) => key(&[b"gi", &gi.to_le_bytes()]),
        ParsedSeqId::Pdb(pdb) => key(&[b"pdb", &pdb.mol]),
        _ => {
            let text = id
                .textseq()
                .expect("the other alternatives are text Seq-ids");
            if let Some(accession) = &text.accession {
                key(&[b"accession", &accession.to_ascii_uppercase()]);
            }
            if let Some(name) = &text.name {
                key(&[b"name", &name.to_ascii_uppercase()]);
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Patent_seq_id.cpp:53-56
// ```c++
// bool CPatent_seq_id::Match(const CPatent_seq_id& psip2) const
// {
//     return GetSeqid() == psip2.GetSeqid()  &&  GetCit().Match(psip2.GetCit());
// }
// ```
// NCBI reference (598d8ae6): c++/src/objects/biblio/Id_pat.cpp:55-79
// ```c++
// bool CId_pat::Id_Match(const C_Id& id1, const C_Id& id2)
// {
// 	C_Id::E_Choice type1 = id1.Which();
// 	C_Id::E_Choice type2 = id2.Which();
//
// 	if (type1 != type2)
// 		return false;
//
//     switch ( type1 ) {
//     case C_Id::e_Number:
//         return AStrEquiv(id1.GetNumber(), id2.GetNumber(), PNocase());
//     case C_Id::e_App_number:
//         return AStrEquiv(id1.GetApp_number(), id2.GetApp_number(), PNocase());
//     default:
//         return false;
//     }
// }
//
//
// // comparison function
// bool CId_pat::Match(const CId_pat& idp2) const
// {
//     return AStrEquiv(GetCountry(), idp2.GetCountry(), PNocase()) &&
//         Id_Match(GetId(), idp2.GetId());
// }
// ```
// `doc-type` takes no part.
fn patent_match(a: &PatentSeqId, b: &PatentSeqId) -> bool {
    a.seqid == b.seqid
        && a.cit.country.eq_ignore_ascii_case(&b.cit.country)
        && match (&a.cit.id, &b.cit.id) {
            (PatentNumber::Number(x), PatentNumber::Number(y))
            | (PatentNumber::AppNumber(x), PatentNumber::AppNumber(y)) => x.eq_ignore_ascii_case(y),
            _ => false,
        }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Textseq_id.cpp:143-174
// ```c++
// bool CTextseq_id::Match(const CTextseq_id& tsip2) const
// {
//     // Check Accessions first
//     if (IsSetAccession()  &&  tsip2.IsSetAccession()) {
//         if ( PNocase().Equals(GetAccession(), tsip2.GetAccession()) ) {
//             if (IsSetVersion()  &&  tsip2.IsSetVersion()) {
//                 return GetVersion() == tsip2.GetVersion();
//             } else {
//                 return true;
//             }
//         } else {
//             return false;
//         }
//     }
//
//     // then try name
//     if (IsSetName()  &&  tsip2.IsSetName()) {
//         if ( PNocase().Equals(GetName(), tsip2.GetName()) ) {
//             if (IsSetVersion()  &&  tsip2.IsSetVersion()) {
//                 return GetVersion() == tsip2.GetVersion();
//             }
//             else {
//                 return true;
//             }
//         } else {
//             return false;
//         }
//     }
//
//     //nothing to compare
//     return false;
// }
// ```
fn textseq_match(a: &TextseqId, b: &TextseqId) -> bool {
    let versions = || match (a.version, b.version) {
        (Some(x), Some(y)) => x == y,
        _ => true,
    };
    if let (Some(x), Some(y)) = (&a.accession, &b.accession) {
        return x.eq_ignore_ascii_case(y) && versions();
    }
    if let (Some(x), Some(y)) = (&a.name, &b.name) {
        return x.eq_ignore_ascii_case(y) && versions();
    }
    false
}

// NCBI reference (598d8ae6): c++/src/objects/general/Dbtag.cpp:158-163
// ```c++
// bool CDbtag::Match(const CDbtag& dbt2) const
// {
//     if (! PNocase().Equals(GetDb(), dbt2.GetDb()))
//         return false;
//     return ((GetTag()).Match((dbt2.GetTag())));
// }
// ```
fn dbtag_match(a: &Dbtag, b: &Dbtag) -> bool {
    a.db.eq_ignore_ascii_case(&b.db) && object_id_match(&a.tag, &b.tag)
}

// NCBI reference (598d8ae6): c++/src/objects/general/Object_id.cpp:93-161
// ```c++
// CObject_id::E_Choice CObject_id::GetIdType(TId8& value) const
// {
//     switch ( Tparent::Which() ) {
//     case e_Id:
//         value = Tparent::GetId();
//         return e_Id;
//     case e_Str:
//     {
//         // try to treat the string as a valid integer the same way as SetStrOrId() does
//         // but also allowing any Int8 values
//         const string& str = GetStr();
//         value = NStr::StringToInt8(str, NStr::fConvErr_NoThrow|NStr::fConvErr_NoErrMessage);
//         if ( !value ) {
//             if ( errno ) {
//                 // not convertible to integer
//                 return CObject_id::e_Str;
//             }
//             // converted to 0
//             if ( str.size() != 1 ) {
//                 // leading zeroes are not allowed
//                 return CObject_id::e_Str;
//             }
//             // valid zero as a string
//             return CObject_id::e_Id;
//         }
//         if ( value > 0 ) {
//             // non-zero positive value
//             if ( str[0] == '0' || str[0] == '+' ) {
//                 // redundant '+' or leading zeroes are not allowed
//                 value = 0;
//                 return CObject_id::e_Str;
//             }
//         }
//         else {
//             // non-zero negative value
//             if ( str[1] == '0' ) {
//                 // leading zeroes are not allowed
//                 value = 0;
//                 return CObject_id::e_Str;
//             }
//         }
//         // valid non-zero value as a string
//         return CObject_id::e_Id;
//     }
//     ...
// // match for identity
// int CObject_id::Compare(const CObject_id& oid2) const
// {
//     TId8 value, value2;
//     E_Choice type = GetIdType(value);
//     E_Choice type2 = oid2.GetIdType(value2);
//     if ( int diff = type - type2 ) {
//         return diff;
//     }
//     switch ( type ) {
//     case e_Id:
//         return value < value2? -1: (value > value2);
//     case e_Str:
//         return PNocase().Compare(GetStr(), oid2.GetStr());
// ```
fn object_id_match(a: &ObjectId, b: &ObjectId) -> bool {
    match (object_id_number(a), object_id_number(b)) {
        (Some(x), Some(y)) => x == y,
        (None, None) => match (a, b) {
            (ObjectId::Str(x), ObjectId::Str(y)) => x.eq_ignore_ascii_case(y),
            _ => false,
        },
        _ => false,
    }
}

/// The value of an Object-id that `GetIdType` reads as an integer: an `id`, or a `str` that
/// is a whole decimal `Int8` without a `+` or leading zeros (`NStr::StringToInt8` with the
/// default flags takes no spaces).
fn object_id_number(tag: &ObjectId) -> Option<i64> {
    match tag {
        ObjectId::Id(id) => Some(i64::from(*id)),
        ObjectId::Str(text) => {
            let text = std::str::from_utf8(text).ok()?;
            let digits = text.strip_prefix('-').unwrap_or(text);
            if digits.is_empty()
                || !digits.bytes().all(|c| c.is_ascii_digit())
                || (digits.len() > 1 && digits.starts_with('0'))
                || (digits == "0" && text.len() != 1)
            {
                return None;
            }
            text.parse::<i64>().ok()
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/PDB_seq_id.cpp:53-67
// ```c++
// bool CPDB_seq_id::Match(const CPDB_seq_id& psip2) const
// {
//     if (IsSetChain() && psip2.IsSetChain()) {
//         if (GetChain() != psip2.GetChain()) {
//             return false;
//         }
//     }
//     if (IsSetChain_id() && psip2.IsSetChain_id()) {
//         if ( PCase().Compare(GetChain_id(), psip2.GetChain_id())) {
//             return false;
//         }
//     }
//     return
//         PCase().Equals(GetMol(), psip2.GetMol());
// }
// ```
fn pdb_match(a: &PdbSeqId, b: &PdbSeqId) -> bool {
    if let (Some(x), Some(y)) = (a.chain, b.chain) {
        if x != y {
            return false;
        }
    }
    if let (Some(x), Some(y)) = (&a.chain_id, &b.chain_id) {
        if x != y {
            return false;
        }
    }
    a.mol == b.mol
}
