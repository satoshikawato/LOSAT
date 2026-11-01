//! A database sequence with parsed Seq-ids in the reports: what each format prints of its
//! definition lines (`seqid::DeflineSet`), chosen and printed by the functions of `seqid`.
//!
//! NCBI's formatters find a subject again by its Seq-align id: the scope gives the Bioseq of
//! the definition line that holds the id (`DeflineSet::bioseq_defline`, the "target") with
//! the whole set as a user object. The ids and title of the target are the Bioseq's (tabular
//! `sseqid`, `sacc`, `saccver`, the description titles of outfmt 0, the XML `Hit_accession`);
//! the formats that show every definition line read the set (the outfmt 0 alignment heading,
//! the XML `Hit_def`, the XML2 and JSON `HitDescr` list, `stitle`).
//!
//! A database without parsed Seq-ids keeps its own path (the title and `gnl|BL_ORD_ID|n`,
//! `algorithm::common::database::read_title`); `SeqDb::read_header` tells the two apart.
//!
//! LOSAT writes the names of the sequence it searched. NCBI's formatters look the Seq-align id
//! up again, and where another sequence of the database has an id that the lookup also matches
//! (the same accession in two sequences, PLAN §9.6 of line DBID) they can write that sequence's
//! names; the test databases leave such records out (`tests/db_search/ids/make_ids.py`), and
//! the check on NCBI's distributed databases is ID4's. The two-letter PDB chains of one
//! molecule (`A`, `AA`, `B`) are not such a case (probe of the ID2a review).

use crate::seqid::defline_set::{
    alignment_line_id, compare, description_id, seq_id_list, tabular_id, Sic, TabularId,
};
use crate::seqid::model::{ObjectId, PatentNumber};
use crate::seqid::rank::{find_best, text_score};
use crate::seqid::{label, BlastDefline, DeflineSet, ParsedSeqId};

use super::blast_xml::HitTaxonomy;
use super::defline::{
    bioseq_title_from_source, generate_bioseq_defline, BioseqIds, Defline, UnknownEncoding,
};

/// A reported database sequence with parsed Seq-ids: its definition lines, its Seq-align id
/// and the definition line of its Bioseq.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParsedSubject {
    set: DeflineSet,
    alignment_id: ParsedSeqId,
    /// `DeflineSet::bioseq_defline` of `alignment_id`.
    target: usize,
    /// Whether NCBI's formatters cannot find the sequence by its Seq-align id (`fetched`).
    unfetchable: bool,
    /// `-show_gis`: the outfmt 0 description line and heading show it (`description`, `heading`).
    show_gis: bool,
    /// The strings of the Seq-aligns' `use_this_seqid` (`under_positive_list`); empty without a
    /// positive list.
    use_this_seqid: Vec<Vec<u8>>,
    /// `tabular_saccver`, kept for the reports of later query batches.
    saccver: std::sync::OnceLock<std::sync::Arc<[u8]>>,
}

/// The `use_this_seqid` strings of a database subject, `subject` when it has parsed Seq-ids: a
/// sequence without them has a `general` id of `BL_ORD_ID`, which is neither a text Seq-id nor a
/// GI (`ParsedSubject::under_positive_list`).
pub fn use_this_seqid(subject: Option<&ParsedSubject>) -> &[Vec<u8>] {
    subject.map_or(&[], ParsedSubject::use_this_seqid)
}

/// `AdjustSeqIdType` with `eGetId_ForceAcc` (`ParsedSubject::sam_id`) of a Bioseq with the ids
/// `ids` and the Seq-align id `alignment_id`.
pub fn forced_accession<'a>(
    alignment_id: &'a ParsedSeqId,
    ids: &'a [ParsedSeqId],
) -> &'a ParsedSeqId {
    let acc_ver =
        |id: &crate::seqid::model::TextseqId| id.accession.is_some() && id.version.is_some();
    if alignment_id.textseq().is_some_and(acc_ver) {
        return alignment_id;
    }
    let mut found = None;
    for id in ids {
        if let Some(text) = id.textseq() {
            found = Some(id);
            if acc_ver(text) {
                break;
            }
        }
    }
    found.unwrap_or(alignment_id)
}

/// One `HitDescr` of the XML2 and JSON formats (`s_SeqAlignToXMLHit`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HitDescription {
    pub id: Vec<u8>,
    /// Not written when empty.
    pub accession: Vec<u8>,
    /// Not written when the definition line has no title member.
    pub title: Option<Vec<u8>>,
    /// Not written when 0.
    pub taxid: i32,
    /// What the taxonomy database said about `taxid` (`Unnamed` where nothing was asked); a
    /// `Named` scientific name is written after the taxid.
    pub taxonomy: HitTaxonomy,
}

impl ParsedSubject {
    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:406-415
    // ```c++
    //     int oid = -1;
    //     ...
    //     if (m_BlastDb->SeqidToOid(*seq_id, oid)) {
    // ```
    // NCBI's formatters find the subject's sequence again from its Seq-align id. For a gibbsq
    // id, a `general` id of the database `PIG`, a named-annotation-track id and a `general`
    // id `BL_ORD_ID` that makeblastdb stored as written, the lookup finds no sequence or
    // another one: NCBI 2.17.0 stops with a CCoreException in outfmt 6 (exit status 255) and
    // writes the title `Unknown`, or another sequence's, in outfmt 0 and 5 (probe of line
    // DBID, DECISIONS 8 (3); `LOSAT/tests/seqid/fasta/variants.fa`). In a version 5 database
    // (LMDB) a local integer id (`lcl|123`) is not found either (the same failures; probe of
    // ID2a, `blastdb-ids-20261021/R/probe`), where version 4 (ISAM) finds it. LOSAT refuses
    // these subjects instead of repeating a defect of NCBI's. So it does where the definition
    // line of the Bioseq is found by the content of a local id alone (`lcl|777` before the
    // definition line of the Seq-align id `gi|777`, `DeflineSet::bioseq_defline`): NCBI 2.17.0
    // finds no Bioseq then, in version 4 and 5 (outfmt 6 stops with a CCoreException, XML2 writes
    // the Seq-align id alone; probe of the ID2a review).
    /// The subject of the definition lines `set` of a database of format version 5
    /// (`version5`) or 4. `Err` names what LOSAT does not write as NCBI: a stored title or
    /// Seq-id with a byte outside ' ' ..= '~' (NCBI writes an error message each time it reads
    /// the definition line) and a definition line without a Seq-id (NCBI's formatters take the
    /// best of none; makeblastdb writes none such).
    pub fn new(set: DeflineSet, version5: bool) -> Result<Self, String> {
        if set
            .deflines
            .iter()
            .any(|defline| defline.replaced_title_bytes > 0 || defline.replaced_id_bytes > 0)
        {
            return Err("has a stored title or Seq-id with a byte outside ' ' to '~' (NCBI BLAST+ shows it as '#' and writes an error message each time it reads the definition line)".into());
        }
        if set.deflines.iter().any(|defline| defline.ids.is_empty()) {
            return Err("has a definition line without a Seq-id".into());
        }
        let Some(alignment_id) = set.alignment_id().cloned() else {
            return Err("has no Seq-id".into());
        };
        let Some(target) = set.bioseq_defline(&alignment_id) else {
            return Err("has no definition line that holds its Seq-align id".into());
        };
        let unfetchable = unfetchable(&alignment_id, version5)
            || set
                .display_defline(&alignment_id)
                .is_some_and(|shown| shown.index != target);
        Ok(Self {
            set,
            alignment_id,
            target,
            unfetchable,
            show_gis: false,
            use_this_seqid: Vec::new(),
            saccver: std::sync::OnceLock::new(),
        })
    }

    /// The subject as the outfmt 0 report of a run with or without `-show_gis` shows it.
    pub fn showing_gis(mut self, show_gis: bool) -> Self {
        self.show_gis = show_gis;
        self
    }

    /// The subject as a database with a positive list (`positive_list`,
    /// `SeqDb::has_positive_list`) reports it: its Seq-aligns name the ids of the definition
    /// lines the list keeps. `use_gi` is whether the Seq-align id is a GI.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqinfosrc_aux.cpp:113-150
    /// ```c++
    /// void GetFilteredRedundantSeqids(const IBlastSeqInfoSrc & seqinfo_src,
    ///                              int                      oid,
    ///                              vector<string>            & seqids,
    /// 			     bool			use_gi)
    /// {
    ///     seqids.resize(0);
    ///
    ///     if (! seqinfo_src.HasGiList()) {
    ///         return;
    ///     }
    ///
    ///     list< CRef<CSeq_id> > seqid_list = seqinfo_src.GetId(oid);
    ///
    ///     ITERATE(list< CRef<CSeq_id> >, id, seqid_list) {
    /// 	if (use_gi)
    /// 	{
    ///         	if ((**id).IsGi())
    /// 		{
    /// 			string sid_string = "gi:" + (*id)->GetSeqIdString(true);
    /// 			seqids.push_back(sid_string);
    /// 		}
    /// 	}
    /// 	else
    /// 	{
    /// 		const CTextseq_id* tsid = (*id)->GetTextseq_Id();
    /// 		if (tsid && tsid->IsSetAccession())
    /// 		{
    /// 			string sid_string = "seqid:" + (*id)->GetSeqIdString(true);
    /// 			seqids.push_back(sid_string);
    /// 		}
    /// 		else if ((*id)->IsPdb())
    /// 		{
    /// 			string sid_string = "seqid:" + (*id)->GetSeqIdString(true);
    /// 			seqids.push_back(sid_string);
    /// 		}
    /// 	}
    ///     }
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1604-1605
    /// ```c++
    ///         // Get SeqIds for entrez query restriction.
    /// 	vector<string> seqid_list;
    ///         GetFilteredRedundantSeqids(*seqinfo_src, hsp_list->oid, seqid_list, subject_id->IsGi());
    /// ```
    pub fn under_positive_list(mut self, positive_list: bool) -> Self {
        if !positive_list {
            return self;
        }
        let use_gi = self.alignment_id.gi().is_some();
        for id in self.set.ids() {
            let prefix: &[u8] = if use_gi {
                if id.gi().is_none() {
                    continue;
                }
                b"gi:"
            } else {
                if !(id.textseq().is_some_and(|text| text.accession.is_some())
                    || matches!(id, ParsedSeqId::Pdb(_)))
                {
                    continue;
                }
                b"seqid:"
            };
            let mut string = prefix.to_vec();
            crate::seqid::label::content(id, true, &mut string);
            self.use_this_seqid.push(string);
        }
        self
    }

    /// The strings of the Seq-aligns' `use_this_seqid` (`under_positive_list`).
    pub fn use_this_seqid(&self) -> &[Vec<u8>] {
        &self.use_this_seqid
    }

    /// `Err` for a subject whose sequence NCBI's formatters cannot find by its Seq-align id
    /// (see `new`): the formats that read the subject's Bioseq (outfmt 0, 5, 6, 7, 13 to 16 and
    /// SAM, 17) refuse it; the Seq-annots and the archive (8, 9, 11, 12) write the Seq-align id
    /// alone.
    pub fn fetched(&self) -> Result<(), String> {
        if !self.unfetchable {
            return Ok(());
        }
        Err(format!(
            "has the Seq-align id {}, by which NCBI BLAST+'s formatters do not find the sequence (NCBI BLAST+ 2.17.0 writes Unknown or another sequence's names, or stops)",
            String::from_utf8_lossy(&self.alignment_fasta())
        ))
    }

    /// The Seq-align id (outfmt 8, 9, 11, 12; `DeflineSet::alignment_id`).
    pub fn alignment_id(&self) -> &ParsedSeqId {
        &self.alignment_id
    }

    // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:131-141
    // ```c++
    // void CSAM_CIGAR_Formatter::AdjustSeqIdType(CConstRef<CSeq_id>& id)
    // {
    //     ...
    //     CSeq_id_Handle forced_id = sequence::GetId(*id, *scope, force_type);
    //     if (forced_id) {
    //         id.Reset(forced_id.GetSeqId());
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:2743-2757 (`GetAccVer`, which
    // `GetId(eGetId_ForceAcc)` calls: the id itself where it is an accession with a version,
    // else `x_GetAccVer` of the Bioseq's ids)
    // ```c++
    //     if ( !(flags & CScope::fForceLoad) ) {
    //         if ( idh.IsAccVer() ) {
    //             return idh;
    //         }
    //     }
    //     ...
    //             CSeq_id_Handle ret = CScope::x_GetAccVer(info->GetIds());
    // ```
    // NCBI reference (598d8ae6): c++/src/objmgr/scope.cpp:425-445
    // ```c++
    //     ITERATE ( TIds, iter, ids ) {
    //         if ( iter->IsGi() ) {
    //             // optimization to avoid calling GetSeqId() (it may be expensive)
    //             continue;
    //         }
    //         if ( auto seq_id = iter->GetSeqId() ) { // need to hold CConstRef<>
    //             if ( auto text_id = seq_id->GetTextseq_Id() ) {
    //                 // found some text id, it can be returned
    //                 ret = *iter;
    //                 if ( text_id->IsSetAccession() && text_id->IsSetVersion() ) {
    //                     // found acc.ver, stop searching
    //                     break;
    //                 }
    // ```
    /// The Seq-id that NCBI's SAM formatter (outfmt 17) names the subject by
    /// (`AdjustSeqIdType` with `eGetId_ForceAcc`): the Seq-align id where it is an accession
    /// with a version, else the first such id of the subject Bioseq (its gis passed over), else
    /// its last Textseq-id, else the Seq-align id.
    pub fn sam_id(&self) -> &ParsedSeqId {
        forced_accession(&self.alignment_id, &self.target().ids)
    }

    /// The definition lines.
    pub fn deflines(&self) -> &[BlastDefline] {
        &self.set.deflines
    }

    /// The definition line of the subject Bioseq.
    fn target(&self) -> &BlastDefline {
        &self.set.deflines[self.target]
    }

    /// What `GenerateDefline` reads of the subject Bioseq's Seq-ids.
    fn bioseq_ids(&self) -> BioseqIds {
        bioseq_ids(&self.target().ids)
    }

    /// `GenerateDefline` of the subject Bioseq: for the alignment heading
    /// (`leave_prefix_suffix` false) or the description table (true).
    pub fn bioseq_title(
        &self,
        protein: bool,
        leave_prefix_suffix: bool,
    ) -> Result<Defline, UnknownEncoding> {
        let title = self.target().title.as_deref().unwrap_or_default();
        generate_bioseq_defline(title, protein, leave_prefix_suffix, &self.bioseq_ids())
    }

    /// Whether the subject Bioseq has a source: some definition line has a positive taxid
    /// (`CSeqDBVol::x_GetTaxonomy` makes a source of each, seqdbvol.cpp:1136).
    fn has_source(&self) -> bool {
        !self.volume_taxids().is_empty()
    }

    /// The positive taxids of the definition lines in the order in which the volume looks them
    /// up when it fetches the subject Bioseq (`CSeqDBVol::x_GetTaxonomy`; `x_GetTaxDefline`
    /// puts the target first), each once.
    pub fn volume_taxids(&self) -> Vec<i32> {
        let mut taxids = Vec::new();
        let order = std::iter::once(self.target()).chain(
            (0..self.set.deflines.len())
                .filter(|&i| i != self.target)
                .map(|i| &self.set.deflines[i]),
        );
        for defline in order {
            match defline.taxid {
                Some(taxid) if taxid > 0 && !taxids.contains(&taxid) => taxids.push(taxid),
                _ => {}
            }
        }
        taxids
    }

    /// Whether NCBI makes the subject Bioseq's title from its source (`bioseq_title_from_source`)
    /// and that source has a taxid: the title then depends on the taxonomy database.
    pub fn title_from_taxonomy(&self) -> bool {
        let title = self.target().title.as_deref().unwrap_or_default();
        bioseq_title_from_source(title, &self.bioseq_ids()) && self.has_source()
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3964-3972
    // ```c++
    //     // adjust protein partial/organelle/taxname suffix, if necessary
    //     if ( m_IsAA && (! m_MainTitle.empty()) && (! (flags & fLeavePrefixSuffix)) ) {
    //         if (m_Idx) {
    //             x_AdjustProteinTitleSuffixIdx (bsh);
    //         } else {
    //             x_AdjustProteinTitleSuffix (bsh);
    //         }
    //     }
    // ```
    // `x_AdjustProteinTitleSuffix` (:3263-3330) returns at once for a Bioseq without a source.
    // Otherwise a title that ends in ` [<taxname>]` (the source's name from the taxonomy
    // database, empty without it) loses that suffix, a `(organelle)` and a `, partial` before
    // it, and gets the name back as the taxonomy database writes it: `alpha, partial [Homo
    // sapiens]` is written `alpha [Homo sapiens]` (probe of the ID2a review).
    /// Whether NCBI's `GenerateDefline` without `fLeavePrefixSuffix` may rewrite the end of the
    /// subject Bioseq's title (`x_AdjustProteinTitleSuffix`, which LOSAT does not port): a
    /// protein whose Bioseq has a source and whose title ends in `]`.
    pub fn title_suffix_adjusted(&self, protein: bool) -> bool {
        let title = self.target().title.as_deref().unwrap_or_default();
        let end = title
            .iter()
            .rposition(|byte| !b".,;~ ".contains(byte) && !byte.is_ascii_whitespace());
        protein && self.has_source() && end.is_some_and(|end| title[end] == b']')
    }

    /// Whether the outfmt 0 alignment heading shows the subject Bioseq's generated title: some
    /// definition line has no title (`heading`).
    pub fn heading_generates_title(&self) -> bool {
        self.set
            .deflines
            .iter()
            .any(|defline| defline.title.as_deref().is_none_or(<[u8]>::is_empty))
    }

    /// The first word of the subject Bioseq's generated title, which `seq_id_list` puts in the
    /// place of a `gnl|BL_ORD_ID` or `lcl|Subject_` id of the Bioseq; empty when no such id is
    /// there (the title is then not made). `Err` where NCBI's `x_CleanAndCompress` reads past the
    /// end of the title (NCBI crashes) and where NCBI may rewrite its end
    /// (`title_suffix_adjusted`).
    fn title_word(&self, protein: bool) -> Result<Vec<u8>, &'static str> {
        let ids = &self.target().ids;
        if matches!(seq_id_list(ids, b"x"), std::borrow::Cow::Borrowed(_)) {
            return Ok(Vec::new());
        }
        self.generated_title_word(protein)
    }

    /// The first word of the subject Bioseq's generated title (`title_word` without the check
    /// of the target's ids).
    fn generated_title_word(&self, protein: bool) -> Result<Vec<u8>, &'static str> {
        if self.title_suffix_adjusted(protein) {
            return Err("a title ending in ']' whose first word NCBI BLAST+ takes after it rewrites the title's end from the taxids (x_AdjustProteinTitleSuffix)");
        }
        let defline = self
            .bioseq_title(protein, false)
            .map_err(|UnknownEncoding| "a title in no encoding that NCBI BLAST+ recognises")?;
        // A nucleotide title stops at the end of the string (approved exception 2 of
        // PD-LOSAT-NCBI-DEFECTS for BLASTN, TBLASTN and TBLASTX).
        if protein && defline.reads_past_end {
            return Err("a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes)");
        }
        let end = defline
            .text
            .iter()
            .position(|&byte| byte == b' ')
            .unwrap_or(defline.text.len());
        Ok(defline.text[..end].to_vec())
    }

    /// The tabular `sseqid`, `sacc` and `saccver` (`tabular_id` of the subject Bioseq's ids after
    /// `seq_id_list`).
    pub fn tabular_ids(&self, protein: bool) -> Result<[Vec<u8>; 3], &'static str> {
        let word = self.title_word(protein)?;
        let ids = seq_id_list(&self.target().ids, &word);
        let field = |kind| {
            let mut out = Vec::new();
            tabular_id(&ids, kind, &mut out);
            out
        };
        Ok([
            field(TabularId::Full),
            field(TabularId::Accession),
            field(TabularId::AccessionVersion),
        ])
    }

    /// The `saccver` of `tabular_ids`, made once for a subject (`protein` is its database's).
    pub fn tabular_saccver(&self, protein: bool) -> Result<std::sync::Arc<[u8]>, &'static str> {
        if let Some(saccver) = self.saccver.get() {
            return Ok(std::sync::Arc::clone(saccver));
        }
        let [_, _, saccver] = self.tabular_ids(protein)?;
        Ok(std::sync::Arc::clone(
            self.saccver.get_or_init(|| std::sync::Arc::from(saccver)),
        ))
    }

    /// `Err` as `tabular_ids`', which it does not make.
    pub fn check_tabular_ids(&self, protein: bool) -> Result<(), &'static str> {
        self.title_word(protein).map(drop)
    }

    /// The tabular `sgi`: the first GI of the subject Bioseq's ids (`FindGi`; `seq_id_list` only
    /// replaces ids that are no GI), 0 for none.
    pub fn tabular_gi(&self) -> i64 {
        self.target().gi()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:544-556
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
    /// The tabular ids of every definition line (`x_SetSubjectIds`: `sallseqid`, `sallgi`,
    /// `sallacc`), in the stored order: each line's ids after `seq_id_list` with the first word
    /// of the subject Bioseq's title, as `tabular_ids` makes them. `Err` as `tabular_ids`'.
    pub fn all_tabular_ids(
        &self,
        protein: bool,
    ) -> Result<Vec<super::tabular::DeflineIds>, &'static str> {
        // `GetSeqIdList` replaces the ids of each definition line, so the word is made when any
        // line has an id to replace, not only the target's.
        let replaced =
            self.set.deflines.iter().any(|defline| {
                matches!(seq_id_list(&defline.ids, b"x"), std::borrow::Cow::Owned(_))
            });
        let word = if replaced {
            self.generated_title_word(protein)?
        } else {
            Vec::new()
        };
        Ok(self
            .set
            .deflines
            .iter()
            .map(|defline| {
                let ids = seq_id_list(&defline.ids, &word);
                let field = |kind| {
                    let mut out = Vec::new();
                    tabular_id(&ids, kind, &mut out);
                    out
                };
                super::tabular::DeflineIds {
                    id: field(TabularId::Full),
                    gi: defline.gi(),
                    accession: field(TabularId::Accession),
                }
            })
            .collect())
    }

    /// The title of the first definition line, empty for none (`stitle`; the writers show an
    /// empty one as `N/A`).
    pub fn first_title(&self) -> &[u8] {
        self.set
            .deflines
            .first()
            .and_then(|defline| defline.title.as_deref())
            .unwrap_or_default()
    }

    /// The start of the outfmt 0 description line (`description_id` of `GetDisplayIds`' id) and
    /// its title: `GenerateDefline` of the subject Bioseq with `fLeavePrefixSuffix` (`Unknown`
    /// where NCBI's `HtmlDecode` throws), after two spaces where `gi|N|` was written.
    pub fn description(&self, protein: bool) -> (Vec<u8>, Defline) {
        let (start, wrote_gi) = self.description_id();
        let mut title = self.bioseq_title(protein, true).unwrap_or(Defline {
            text: b"Unknown".to_vec(),
            reads_past_end: false,
        });
        if wrote_gi {
            title.text.splice(0..0, *b"  ");
        }
        (start, title)
    }

    /// The start of the outfmt 0 description line (`description`).
    pub fn description_start(&self) -> Vec<u8> {
        self.description_id().0
    }

    /// `description`'s start and whether `gi|N|` was written.
    fn description_id(&self) -> (Vec<u8>, bool) {
        let mut start = Vec::new();
        let mut wrote_gi = false;
        if let Some(shown) = self.set.display_defline(&self.alignment_id) {
            let defline = &self.set.deflines[shown.index];
            if let Some(id) = defline.display_id() {
                let gi = if shown.found { defline.gi() } else { 0 };
                wrote_gi = description_id(id, gi, self.show_gis, &mut start);
            }
        }
        (start, wrote_gi)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2380-2384, 2451-2461
    // ```c++
    //                 if(alnDispParams) {
    //                     numBdl++;
    //                     if(isFirst){
    //                         out << ">";
    //                     } else{
    //                         out << " ";
    //     ...
    //                     if (out.tellp() > 1L) {
    //                         out << " ";
    //                     }
    // 					if(!alnDispParams->title.empty()) {
    //                         s_WrapOutputLine(out, (m_AlignOption&eHtml) ?
    //                                          CHTMLHelper::
    //                                          HTMLEncode(alnDispParams->title) :
    //                                          alnDispParams->title);
    //                     }
    //                     out<<"\n";
    //                     isFirst = false;
    // ```
    // `out` holds the heading of this subject only; the title of a definition line is its own,
    // or the subject Bioseq's generated title when it is empty (`x_FillDispParams`).
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3166-3180
    // ```c++
    //     ITERATE(CSeq_align_set::Tdata, iter, alnset.Get()){
    //         const CSeq_align::TSegs& seg = (*iter)->GetSegs();
    //         if(seg.Which() == CSeq_align::C_Segs::e_Std){
    //             if(seg.GetStd().size() > 1){
    //                 //has more than one stdseg. Need to seperate as each
    //                 //is a distinct HSP
    //                 ITERATE (CSeq_align::C_Segs::TStd, iterStdseg, seg.GetStd()){
    //                     CRef<CSeq_align> aln(new CSeq_align);
    //                     if((*iterStdseg)->IsSetScores()){
    //                         aln->SetScore() = (*iterStdseg)->GetScores();
    //                     }
    //                     aln->SetSegs().SetStd().push_back(*iterStdseg);
    //                     alnSetRef->Set().push_back(aln);
    //                 }
    //
    //             } else {
    //                 alnSetRef->Set().push_back(*iter);
    //             }
    // ```
    /// The definition-line lines of the outfmt 0 alignment of the subject (each ending with a
    /// newline; `Length=` follows): the definition lines (`heading_shows`), the first after `>`
    /// and the others after a space, with its id (`alignment_line_id`) and its title wrapped at
    /// 60 bytes. `every_line`: the Seq-align has no `ext` when the report is made, so nothing
    /// limits the lines: a Seq-align of an ungapped search with two or more HSPs, which the
    /// formatter splits into new ones (`PrepareBlastUngappedSeqalign` above).
    pub fn heading(&self, protein: bool, every_line: bool) -> Vec<u8> {
        let mut out = Vec::new();
        let mut generated = None;
        let mut first = true;
        for defline in &self.set.deflines {
            if !every_line && !self.heading_shows(defline) {
                continue;
            }
            let start = out.len();
            out.push(if first { b'>' } else { b' ' });
            alignment_line_id(defline, self.show_gis, &mut out);
            if !first || out.len() - start > 1 {
                out.push(b' ');
            }
            first = false;
            let title = match defline.title.as_deref() {
                Some(title) if !title.is_empty() => title,
                _ => generated.get_or_insert_with(|| {
                    self.bioseq_title(protein, false)
                        .map(|defline| defline.text)
                        .unwrap_or_default()
                }),
            };
            wrap_title(title, &mut out);
            out.push(b'\n');
        }
        out
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2196-2202
    // ```c++
    //     CRef<CSeq_id> wid = FindBestChoice(ids, CSeq_id::WorstRank);
    //     TGi gi_in_use_this_gi = ZERO_GI;
    //     bool isGiList = false;
    //     bool match = CAlignFormatUtil::MatchSeqInSeqList(gi, wid, use_this_seqid,&isGiList);
    //     if(match && isGiList) gi_in_use_this_gi = gi;
    //
    //     if(use_this_seqid.empty() || match) {
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4648-4664
    // ```c++
    // bool CAlignFormatUtil::MatchSeqInSeqList(TGi cur_gi, CRef<CSeq_id> &seqID, list<string> &use_this_seq,bool *isGiList)
    // {
    //     bool found = false;
    //     bool isGi = false;
    //
    //     string curSeqID = CAlignFormatUtil::GetLabel(seqID,true); //uses GetSeqIdString(true)
    //     ITERATE(list<string>, iter_seq, use_this_seq){
    //         isGi = false;
    //         string useThisSeq = s_UseThisSeqToTextSeqID(*iter_seq, isGi);
    //         if((isGi && cur_gi == NStr::StringToNumeric<TGi>((useThisSeq))) || (!isGi && curSeqID == useThisSeq)){
    //             found = true;
    //             break;
    //          }
    //     }
    // ```
    /// Whether the outfmt 0 alignment heading shows `defline`: under a positive list only a line
    /// that `use_this_seqid` names (a line of a local id alone is not shown, nor one without a
    /// GI when the strings are GIs; probe G05 of line DI-8).
    fn heading_shows(&self, defline: &BlastDefline) -> bool {
        if self.use_this_seqid.is_empty() {
            return true;
        }
        // `CAlignFormatUtil::GetLabel`: the tag of a `general` id, else `GetSeqIdString(true)`.
        let mut label_of_id = Vec::new();
        match defline.display_id() {
            Some(id @ ParsedSeqId::General(_)) => {
                label::content_general_db_is_content(id, &mut label_of_id);
            }
            Some(id) => label::content(id, true, &mut label_of_id),
            None => {}
        }
        let gi = defline.gi().to_string();
        self.use_this_seqid
            .iter()
            .any(|string| match string.strip_prefix(b"gi:") {
                Some(value) => value == gi.as_bytes(),
                None => string.strip_prefix(b"seqid:") == Some(label_of_id.as_slice()),
            })
    }

    /// The outfmt 5 `Hit_id` and `Hit_def` (`DeflineSet::hit_id_and_def` with the GIs shown;
    /// an empty `Hit_def` is written `No definition line`) and `Hit_accession` (the label of
    /// the `Score`-best id of the subject Bioseq, without its version; a `general` id's tag).
    pub fn xml_hit(&self) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
        let (id, def) = self.set.hit_id_and_def(true);
        let mut accession = Vec::new();
        if let Some(best) = find_best(&self.target().ids, text_score) {
            label::format_label(best, false, &mut accession);
        }
        (id, def, accession)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:365-388
    // ```c++
    //         	ITERATE(list<CRef<CBlast_def_line> >, itr, bdlRef->Get()) {
    //         		const CBlast_def_line & defline = **itr;
    //         		CRef<blastxml2::CHitDescr> hit_exp(new blastxml2::CHitDescr);
    //         		hit_exp->SetId(CShowBlastDefline::GetSeqIdListString(defline.GetSeqid(), true));
    //
    //         		 CRef<CSeq_id> best_id = FindBestChoice(defline.GetSeqid(), CSeq_id::Score);
    //         		 CSeq_id_Handle id_handle = CSeq_id_Handle::GetHandle(*best_id);
    //         		 string accession = CAlignFormatUtil::GetLabel(id_handle.GetSeqId());
    //         		 if(accession != kEmptyStr)
    //         			 hit_exp->SetAccession(accession);
    //
    //         		if(defline.IsSetTitle())
    //         			hit_exp->SetTitle(defline.GetTitle());
    //
    //         		if(defline.IsSetTaxid() && defline.GetTaxid() != ZERO_TAX_ID) {
    //         			TTaxId tax_id = defline.GetTaxid();
    //        				hit_exp->SetTaxid(tax_id);
    // ```
    /// The `HitDescr` of each definition line (XML2 and JSON), without the scientific names.
    pub fn hit_descriptions(&self) -> Vec<HitDescription> {
        self.set
            .deflines
            .iter()
            .map(|defline| {
                let mut id = Vec::new();
                label::seqid_list_string(&defline.ids, true, &mut id);
                let mut accession = Vec::new();
                if let Some(best) = defline.display_id() {
                    label::format_label(best, false, &mut accession);
                }
                HitDescription {
                    id,
                    accession,
                    title: defline.title.clone(),
                    taxid: defline.taxid.unwrap_or(0),
                    taxonomy: HitTaxonomy::Unnamed,
                }
            })
            .collect()
    }

    /// The FASTA form of the Seq-align id (`CSeq_id::AsFastaString`): the id of the `HitDescr`
    /// that NCBI writes in place of the rest when the taxonomy database does not know a taxid
    /// (blastxml2_format.cpp:422-428), and the id that LOSAT's refusals name.
    pub fn alignment_fasta(&self) -> Vec<u8> {
        let mut out = Vec::new();
        label::fasta(&self.alignment_id, &mut out);
        out
    }
}

/// Whether NCBI's formatters cannot find a sequence by the Seq-align id `id` in a database of
/// version 5 (`version5`) or 4 (`ParsedSubject::new`).
fn unfetchable(id: &ParsedSeqId, version5: bool) -> bool {
    match id {
        ParsedSeqId::Gibbsq(_) | ParsedSeqId::NamedAnnotTrack(_) => true,
        ParsedSeqId::Local(ObjectId::Id(_)) => version5,
        ParsedSeqId::General(dbtag) => dbtag.db == b"PIG" || dbtag.db == b"BL_ORD_ID",
        _ => false,
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/accguide2.txt (the rows of the RefSeq
// kinds that `x_SetFlags` tests; `3+6` is a prefix of three characters and six digits)
// ```
// 3+6  AC_      refseq_chromosome
// 3+6  NC_      refseq_chromosome
// 3+6  NM_      refseq_mrna
// 3+9  NM_      refseq_mrna
// 3+6  NR_      refseq_ncrna
// 3+6  XM_      refseq_mrna_predicted
// 3+9  XM_      refseq_mrna_predicted
// ```
/// Whether `CSeq_id::IdentifyAccession` gives `accession` a RefSeq chromosome, mRNA or ncRNA
/// type (`BioseqIds::refseq_source_title`).
fn refseq_source_accession(accession: &[u8]) -> bool {
    if accession.len() < 3 || accession[2] != b'_' {
        return false;
    }
    let digits = &accession[3..];
    if !digits.iter().all(u8::is_ascii_digit) {
        return false;
    }
    let prefix = [
        accession[0].to_ascii_uppercase(),
        accession[1].to_ascii_uppercase(),
    ];
    match &prefix {
        b"AC" | b"NC" | b"NR" => digits.len() == 6,
        b"NM" | b"XM" => digits.len() == 6 || digits.len() == 9,
        _ => false,
    }
}

/// `BioseqIds` of the Seq-ids `ids` of a Bioseq.
pub(crate) fn bioseq_ids(ids: &[ParsedSeqId]) -> BioseqIds {
    let mut flags = BioseqIds::default();
    for id in ids {
        match id {
            ParsedSeqId::Tpg(text) | ParsedSeqId::Tpe(text) | ParsedSeqId::Tpd(text) => {
                flags.third_party = true;
                flags.refseq_source_title |= text
                    .accession
                    .as_deref()
                    .is_some_and(refseq_source_accession);
            }
            ParsedSeqId::Other(text)
            | ParsedSeqId::Genbank(text)
            | ParsedSeqId::Embl(text)
            | ParsedSeqId::Ddbj(text) => {
                flags.refseq_source_title |= text
                    .accession
                    .as_deref()
                    .is_some_and(refseq_source_accession);
            }
            ParsedSeqId::Pdb(pdb) => {
                let (chain_id, chain) = flags.pdb.get_or_insert_with(|| (Vec::new(), 0));
                if let Some(id) = &pdb.chain_id {
                    *chain_id = id.clone();
                } else if let Some(value) = pdb.chain {
                    *chain = value;
                }
            }
            ParsedSeqId::Patent(patent) => {
                let number = match &patent.cit.id {
                    PatentNumber::Number(number) | PatentNumber::AppNumber(number) => {
                        number.clone()
                    }
                };
                flags.patent = Some((patent.seqid, patent.cit.country.clone(), number));
            }
            _ => {}
        }
    }
    flags
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:336-359
// ```c++
// static void s_WrapOutputLine(CNcbiOstream& out, const string& str)
// {
//     const int line_len = 60;
//     bool do_wrap = false;
//     int length = (int) str.size();
//     if (length > line_len) {
//         for (int i = 0; i < length; i ++){
//             if(i > 0 && i % line_len == 0){
//                 do_wrap = true;
//             }
//             out << str[i];
//             if(do_wrap && isspace((unsigned char) str[i])){
//                 out << "\n";
//                 do_wrap = false;
//             }
//         }
//     } else {
//         out << str;
//     }
// }
// ```
/// `s_WrapOutputLine`: `title`, with a newline after the first white space at or after each
/// multiple of 60 bytes.
pub fn wrap_title(title: &[u8], out: &mut Vec<u8>) {
    let mut do_wrap = false;
    for (i, &c) in title.iter().enumerate() {
        if i > 0 && i % 60 == 0 {
            do_wrap = true;
        }
        out.push(c);
        if do_wrap && (c.is_ascii_whitespace() || c == 11) {
            out.push(b'\n');
            do_wrap = false;
        }
    }
}

/// Whether NCBI's XML and outfmt 0 reports join the hits of `a` and `b` when they follow each
/// other in a hit list: their Seq-align ids match (`CSeq_id::Match`, which is
/// `Compare(..) == e_YES`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:478-487
/// ```c++
///             CSeq_align_set one_subject_alnset;
///             CConstRef<CSeq_id> current_id(s_GetSubjectId(*(*iter)));
///             for ( ; iter != alnset.Get().end(); ++iter) {
///                 CConstRef<CSeq_id> next_id(s_GetSubjectId(*(*iter)));
///                 if (!current_id->Match(*next_id)) {
///                     break;
///                 }
///                 one_subject_alnset.Set().push_back(*iter);
///             }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1750-1757
/// ```c++
///     ITERATE(CSeq_align_set::Tdata, iter, source.Get()){
///         const CSeq_id& cur_id = (*iter)->GetSeq_id(1);
///         if(previous_id.Empty()) {
///             temp =  new CSeq_align_set;
///             temp->Set().push_back(*iter);
///             target.push_back(temp);
///         } else if (cur_id.Match(*previous_id)){
///             temp->Set().push_back(*iter);
/// ```
pub fn joined_in_reports(a: &ParsedSubject, b: &ParsedSubject) -> bool {
    compare(&a.alignment_id, &b.alignment_id) == Sic::Yes
}

// No NCBI counterpart: LOSAT's check that no two sequences of a database are joined
// (`common::database::reports_may_refuse`).
/// Adds to `keys` the keys of the Seq-align id of `subject` (`defline_set::match_keys`): two
/// subjects that the reports join (`joined_in_reports`) share one.
pub fn join_keys(subject: &ParsedSubject, keys: &mut Vec<u64>) {
    crate::seqid::defline_set::match_keys(&subject.alignment_id, keys);
}

/// The keys of the Seq-align id of `subject` as NCBI's scope finds it
/// (`defline_set::scope_keys`).
pub fn scope_keys(subject: &ParsedSubject, keys: &mut Vec<u64>) {
    crate::seqid::defline_set::scope_keys(&subject.alignment_id, keys);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqdb::{SeqDb, SeqType};
    use crate::seqid::model::{Dbtag, PdbSeqId};

    /// The subjects of `tests/db_search/ids/prot_v4` (`make_ids.py`), by the FASTA form of the
    /// first id of their first definition line.
    fn subject(first: &str) -> ParsedSubject {
        let base =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids/prot_v4");
        let db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
        for ordinal in 0..db.num_seqs() {
            let set = db.defline_set(ordinal).unwrap();
            let mut fasta = Vec::new();
            label::fasta(&set.deflines[0].ids[0], &mut fasta);
            if fasta == first.as_bytes() {
                return ParsedSubject::new(set, false).unwrap();
            }
        }
        panic!("no sequence {first}");
    }

    fn text(bytes: &[u8]) -> &str {
        std::str::from_utf8(bytes).unwrap()
    }

    // Approved exception 2 of PD-LOSAT-NCBI-DEFECTS: BLASTN, TBLASTN and TBLASTX stop an outfmt 0
    // title made only of punctuation at the end of the string (NCBI's `x_CleanAndCompress` reads
    // past it and crashes); a nucleotide subject with parsed Seq-ids gets the title of the same
    // stored title without them, and so does the tabular title word. A protein subject is refused.
    #[test]
    fn a_nucleotide_title_of_punctuation_stops_at_its_end() {
        use crate::algorithm::common::database::check_shown_parsed_subject;
        let subject = |id: ParsedSeqId| {
            let defline = BlastDefline {
                title: Some(b", ,".to_vec()),
                ids: vec![id],
                ..BlastDefline::default()
            };
            ParsedSubject::new(
                DeflineSet {
                    deflines: vec![defline],
                },
                true,
            )
            .unwrap()
        };
        let stored = crate::report::defline::generate_defline(b", ,", false, true).unwrap();
        assert!(stored.reads_past_end);
        let accession = subject(ParsedSeqId::Genbank(crate::seqid::model::TextseqId {
            name: None,
            accession: Some(b"AF000050".to_vec()),
            release: None,
            version: Some(1),
        }));
        let (start, title) = accession.description(false);
        assert_eq!(text(&start), "AF000050.1 ");
        assert_eq!(title.text, stored.text);
        check_shown_parsed_subject(&accession, false, false, 0, "BLASTN")
            .expect("a nucleotide title of punctuation is shown");
        assert!(check_shown_parsed_subject(&accession, true, false, 0, "BLASTP").is_err());
        let local = subject(ParsedSeqId::Local(ObjectId::Str(b"Subject_51".to_vec())));
        assert!(local.tabular_ids(true).is_err());
        assert!(local.tabular_saccver(true).is_err());
        assert!(local.check_tabular_ids(true).is_err());
        assert!(local.check_tabular_ids(false).is_ok());
        let [id, _, _] = local.tabular_ids(false).unwrap();
        let word = stored.text.split(|&byte| byte == b' ').next().unwrap();
        assert_eq!(id, word);
    }

    // NCBI BLAST+ 2.17.0 `blastn -db tests/db_search/ids/nucl_v5 -outfmt 17` names the subjects
    // so (`LOSAT/tests/fixtures/db_search/db.id.blastn.17.v5.out`, ordinal order): the Seq-align
    // id `gi|654321` and `gi|2020` (a GI in the second definition line) give the accession of
    // the Bioseq's definition line, `AB000009` has no version, and `lcl|Subject_31` no Textseq-id.
    #[test]
    fn sam_names_the_subject_by_an_accession_of_its_bioseq() {
        let base =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids/nucl_v5");
        let db = SeqDb::open(base.to_str().unwrap(), SeqType::Nucleotide).unwrap();
        let names: Vec<String> = (0..db.num_seqs())
            .map(|ordinal| {
                let subject = ParsedSubject::new(db.defline_set(ordinal).unwrap(), true).unwrap();
                let mut name = Vec::new();
                crate::report::seqid::SeqId::Db(subject.sam_id()).sam_name(&mut name);
                String::from_utf8(name).unwrap()
            })
            .collect();
        assert_eq!(
            names,
            [
                "NM_000001.1",
                "AF000002.3",
                "AJ000003.1",
                "AB000004.1",
                "NC_000005.11",
                "mynuc6",
                "mydb:nuc7",
                "NM_000008.2",
                "AB000009",
                "mynuc_10",
                "BK000012.1",
                "AF000013.1",
                "US1234567_14",
                "AF000015.1",
                "Subject_30",
                "Subject_31",
                "BK000018.1",
                "NC_000019.1",
                "AF000020.1",
                "AF000021.1",
                "AF000022.1",
                "NR_000023.1",
                "ti:123",
                "BN000025.1",
            ]
        );
    }

    fn description(first: &str) -> String {
        let (mut line, title) = subject(first).description(true);
        line.extend_from_slice(&title.text);
        String::from_utf8(line).unwrap()
    }

    // The strings are NCBI BLAST+ 2.17.0's `blastp -db` (outfmt 0 and 6) for these records
    // (`LOSAT/tests/fixtures/db_search/db.id.blastp.0.v4.out`, `db.id.blastp.6c.v4.out`).
    #[test]
    fn titles_that_ncbi_makes_from_the_seq_ids() {
        assert_eq!(
            description("tpg|DAA10001.1|"),
            "DAA10001.1 tpa titled protein"
        );
        assert_eq!(
            description("tpg|DAA10002.1|"),
            "DAA10002.1 unnamed protein product"
        );
        assert_eq!(
            text(&subject("tpg|DAA10002.1|").heading(true, false)),
            ">DAA10002.1 TPA: unnamed protein product\n"
        );
        assert_eq!(description("pdb|2XYZ|B"), "2XYZ_B Chain B");
        assert_eq!(description("pdb|2XYY| "), "2XYY unnamed protein product");
        assert_eq!(
            description("pat|US|7654321|3"),
            "US7654321_3 Sequence 3 from Patent US 7654321"
        );
        assert_eq!(description("ref|NC_100013.1|"), "NC_100013.1 ");
        assert_eq!(
            text(&subject("ref|NC_100013.1|").heading(true, false)),
            ">NC_100013.1 \n"
        );
        assert_eq!(
            description("ref|WP_100014.1|"),
            "WP_100014.1 unnamed protein product"
        );
        assert_eq!(description("lcl|Subject_27"), "subject seven title");
        assert_eq!(
            text(&subject("lcl|Subject_30").heading(true, false)),
            ">Subject_30 unnamed protein product\n"
        );
        assert_eq!(description("pdb|5ABC|D"), "5ABC_D Chain D");
        assert_eq!(
            text(&subject("pdb|5ABC|D").heading(true, false)),
            ">5ABC_D Chain D\n XP_100015.1 second\n"
        );
        assert_eq!(
            text(&subject("lcl|Subject_36").heading(true, false)),
            ">Subject_36 unnamed protein product\n DAA10016.1 tpa second\n"
        );
        // A title from the source: the taxonomy database would name its organism.
        assert!(!subject("ref|NP_100006.1|").title_from_taxonomy());
        assert!(!subject("pdb|2XYZ|B").title_from_taxonomy());
        assert!(subject("gi|1818").title_from_taxonomy());
    }

    #[test]
    fn tabular_ids_replace_subject_ids_by_the_title_word() {
        let ids = |first: &str| {
            subject(first)
                .tabular_ids(true)
                .unwrap()
                .map(|id| String::from_utf8(id).unwrap())
        };
        assert_eq!(ids("lcl|Subject_36"), ["unnamed", "unnamed", "unnamed"]);
        assert_eq!(ids("lcl|Subject_27"), ["subject", "subject", "subject"]);
        assert_eq!(
            ids("lcl|Subject_28"),
            ["tpg|DAA10008.1|", "DAA10008", "DAA10008.1"]
        );
        assert_eq!(ids("lcl|Subject_29"), ["pdb|3ABC|C", "3ABC_C", "3ABC_C"]);
        assert!(subject("lcl|Subject_29").check_tabular_ids(true).is_ok());
        // The `saccver`, made once.
        let once = subject("lcl|Subject_28");
        let first = once.tabular_saccver(true).unwrap();
        assert_eq!(&*first, b"DAA10008.1");
        assert!(std::sync::Arc::ptr_eq(
            &first,
            &once.tabular_saccver(true).unwrap()
        ));
    }

    fn one_defline(ids: Vec<ParsedSeqId>) -> DeflineSet {
        DeflineSet {
            deflines: vec![BlastDefline {
                title: Some(b"t".to_vec()),
                ids,
                ..BlastDefline::default()
            }],
        }
    }

    fn gb(accession: &str, version: i32) -> ParsedSeqId {
        ParsedSeqId::Genbank(crate::seqid::model::TextseqId {
            accession: Some(accession.as_bytes().to_vec()),
            version: Some(version),
            ..Default::default()
        })
    }

    fn lines(lines: Vec<Vec<ParsedSeqId>>) -> ParsedSubject {
        let deflines = lines
            .into_iter()
            .map(|ids| BlastDefline {
                title: Some(b"t".to_vec()),
                ids,
                ..BlastDefline::default()
            })
            .collect();
        ParsedSubject::new(DeflineSet { deflines }, true).unwrap()
    }

    fn strings(subject: &ParsedSubject) -> Vec<String> {
        subject
            .use_this_seqid()
            .iter()
            .map(|string| String::from_utf8(string.clone()).unwrap())
            .collect()
    }

    // Probes G01, G05 of line DI-8 (`GetFilteredRedundantSeqids`, `x_FillAlnDispParams`).
    #[test]
    fn a_positive_list_names_the_text_ids_or_the_gis() {
        let local = || ParsedSeqId::Local(ObjectId::Str(b"loc".to_vec()));
        // Text ids of every line, whichever comes first; a local id is not named.
        let subject = lines(vec![vec![local()], vec![gb("AF000001", 1)]]);
        assert!(subject.use_this_seqid().is_empty());
        let subject = subject.under_positive_list(true);
        assert_eq!(strings(&subject), ["seqid:AF000001.1"]);
        // The heading leaves out the line the strings do not name, and `>` goes to the next.
        assert_eq!(subject.heading(true, false), b">AF000001.1 t\n");
        // An ungapped Seq-align of two HSPs reaches the report without the strings.
        assert_eq!(subject.heading(true, true), b">loc t\n AF000001.1 t\n");
        // A local id alone: no string, every line in the heading.
        let subject = lines(vec![vec![local()]]).under_positive_list(true);
        assert!(subject.use_this_seqid().is_empty());
        assert_eq!(subject.heading(true, false), b">loc t\n");
        // A Seq-align id that is a GI: only GIs are named, and a line without one is left out.
        let subject = lines(vec![
            vec![ParsedSeqId::Gi(888), gb("GG000001", 1)],
            vec![gb("GG000002", 1)],
            vec![ParsedSeqId::Gi(890), gb("GG000003", 1)],
        ])
        .under_positive_list(true);
        assert_eq!(strings(&subject), ["gi:888", "gi:890"]);
        assert_eq!(
            subject.heading(true, false),
            b">GG000001.1 t\n GG000003.1 t\n"
        );
        // A name without an accession is not named; a PDB id is.
        let name_only = ParsedSeqId::Prf(crate::seqid::model::TextseqId {
            name: Some(b"NAME".to_vec()),
            ..Default::default()
        });
        let pdb = ParsedSeqId::Pdb(crate::seqid::model::PdbSeqId {
            mol: b"3DEF".to_vec(),
            chain: None,
            rel: None,
            chain_id: Some(b"B".to_vec()),
        });
        let subject = lines(vec![vec![name_only], vec![pdb]]).under_positive_list(true);
        assert_eq!(strings(&subject), ["seqid:3DEF_B"]);
        // Without a positive list nothing changes.
        let subject =
            lines(vec![vec![local()], vec![gb("AF000001", 1)]]).under_positive_list(false);
        assert!(subject.use_this_seqid().is_empty());
        assert_eq!(subject.heading(true, false), b">loc t\n AF000001.1 t\n");
    }

    #[test]
    fn seq_align_ids_that_ncbi_does_not_find() {
        let local_id = one_defline(vec![ParsedSeqId::Local(ObjectId::Id(123))]);
        assert!(ParsedSubject::new(local_id.clone(), false)
            .unwrap()
            .fetched()
            .is_ok());
        let error = ParsedSubject::new(local_id, true)
            .unwrap()
            .fetched()
            .unwrap_err();
        assert!(error.contains("Seq-align id lcl|123,"), "{error}");
        for id in [
            ParsedSeqId::Gibbsq(123),
            ParsedSeqId::General(Dbtag {
                db: b"PIG".to_vec(),
                tag: ObjectId::Id(123),
            }),
            ParsedSeqId::General(Dbtag {
                db: b"BL_ORD_ID".to_vec(),
                tag: ObjectId::Id(5),
            }),
        ] {
            let subject = ParsedSubject::new(one_defline(vec![id]), false).unwrap();
            assert!(subject.fetched().is_err());
        }
        // A GI in the set makes the GI the Seq-align id.
        let with_gi = one_defline(vec![ParsedSeqId::Gibbsq(123), ParsedSeqId::Gi(7)]);
        assert!(ParsedSubject::new(with_gi, true).unwrap().fetched().is_ok());
        let replaced = DeflineSet {
            deflines: vec![BlastDefline {
                ids: vec![ParsedSeqId::Gi(7)],
                replaced_title_bytes: 1,
                ..BlastDefline::default()
            }],
        };
        assert!(ParsedSubject::new(replaced, false).is_err());
        assert!(ParsedSubject::new(DeflineSet::default(), false).is_err());
    }

    #[test]
    fn pdb_ids_without_a_chain_match_any_chain_of_the_molecule() {
        let pdb = |chain_id: Option<&str>| {
            let id = ParsedSeqId::Pdb(PdbSeqId {
                mol: b"1ABC".to_vec(),
                chain: chain_id.map(|id| i32::from(id.as_bytes()[0])),
                rel: None,
                chain_id: chain_id.map(|id| id.as_bytes().to_vec()),
            });
            ParsedSubject::new(one_defline(vec![id]), false).unwrap()
        };
        assert!(joined_in_reports(&pdb(None), &pdb(Some("AA"))));
        assert!(!joined_in_reports(&pdb(Some("A")), &pdb(Some("B"))));
    }

    #[test]
    fn refseq_accessions_whose_empty_titles_come_from_the_source() {
        for (accession, from_source) in [
            ("NC_000001", true),
            ("AC_000001", true),
            ("NM_000001", true),
            ("NM_000000001", true),
            ("XM_000000001", true),
            ("NR_000001", true),
            ("NR_000000001", false),
            ("NC_0000001", false),
            ("NP_000001", false),
            ("WP_000000001", false),
            ("XR_000001", false),
        ] {
            assert_eq!(
                refseq_source_accession(accession.as_bytes()),
                from_source,
                "{accession}"
            );
        }
    }
}
