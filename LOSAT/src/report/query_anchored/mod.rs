//! The query-anchored alignment section of `-outfmt 1` to `4` (NCBI's `CDisplaySeqalign` with
//! `eMergeAlign`): one query's HSPs merged into one multiple alignment anchored on the query
//! (`aln_mix`, NCBI's `CAlnMix`), read as NCBI's `CAlnVec` reads it (`aln_vec`), and written in
//! lines of `-line_length` columns, each row with its label, its first and last position, and in
//! formats 1 and 2 (`eMasterAnchored`) the letters that it inserts between the query's columns.
//! Formats 1 and 3 (`eShowIdentity`) print `.` for a letter equal to the query's.
//!
//! The rest of a report of these formats is outfmt 0's (`report::pairwise::PairwiseView`): this
//! module writes only what `DisplaySeqalign` writes in place of outfmt 0's alignments.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:633-649 (`s_SetFlags`)
//! ```c++
//!     if (format_type >= CFormattingArgs::eQueryAnchoredIdentities &&
//!         format_type <= CFormattingArgs::eFlatQueryAnchoredNoIdentities) {
//!         flags |= CDisplaySeqalign::eMergeAlign;
//!     }
//!     else {
//!         flags |= CDisplaySeqalign::eShowBlastStyleId |
//!                  CDisplaySeqalign::eShowMiddleLine;
//!     }
//!
//!     if (format_type == CFormattingArgs::eQueryAnchoredIdentities ||
//!         format_type == CFormattingArgs::eFlatQueryAnchoredIdentities) {
//!         flags |= CDisplaySeqalign::eShowIdentity;
//!     }
//!     if (format_type == CFormattingArgs::eQueryAnchoredIdentities ||
//!         format_type == CFormattingArgs::eQueryAnchoredNoIdentities) {
//!         flags |= CDisplaySeqalign::eMasterAnchored;
//!     }
//! ```
//!
//! Line RG T3a writes what BLASTP's Seq-aligns reach: Dense-segs of two rows of width 1 without
//! strands, all in NCBI's first frame group. The frame groups of a translated query
//! (`For reading frame ...`, `s_GetStdsegMasterFrame`), the minus strands, the rows of width 3,
//! Dense-diags and `fForceTranslation` (TBLASTX) are not written here (line RG T3b and T3c): the
//! input (`Segment`) has no strand and no width.

mod aln_mix;
mod aln_vec;

use std::collections::VecDeque;
use std::io::{self, Write};

use aln_vec::AlnMap;

/// A sequence of a query's section, as NCBI's `CAlnMix` tells sequences apart (by their Bioseq
/// handle, `x_IdentifyAlnMixSeq`): the query, or a subject by the caller's key. The caller gives
/// one key per Bioseq: HSPs whose subjects NCBI's scope resolves to one Bioseq share a key.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum SeqKey {
    Query,
    Subject(u32),
}

/// One segment of an HSP's Dense-seg (two rows: the query, then the subject): the query's and
/// the subject's first position (0-based, in the records; `None` for a gap) and its length.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Segment {
    pub(crate) query: Option<u32>,
    pub(crate) subject: Option<u32>,
    pub(crate) len: u32,
}

/// The HSPs of one query's section in the order of its Seq-align set after the sort and the
/// prune of the report (`PruneSeqalign`): each with its subject's key and its segments.
#[derive(Default)]
pub(crate) struct Hsps {
    segments: Vec<Segment>,
    /// Each HSP's subject and the end of its segments in `segments`.
    hsps: Vec<(u32, usize)>,
}

impl Hsps {
    pub(crate) fn is_empty(&self) -> bool {
        self.hsps.is_empty()
    }

    /// Adds an HSP of subject `subject`, its Dense-seg's `starts` (two per segment, the query's
    /// then the subject's, `-1` for a gap, `report::serial::seqalign::Segments::starts`) and
    /// `lens`, each start moved by `query_offset` and `subject_offset` (the remaps of the
    /// ranges, `RemapToQueryLoc` and `s_RemapToSubjectLoc`).
    pub(crate) fn push_dense_seg(
        &mut self,
        subject: u32,
        starts: &[i32],
        lens: &[i32],
        query_offset: u32,
        subject_offset: u32,
    ) {
        let position = |start: i32, offset: u32| (start >= 0).then(|| start as u32 + offset);
        self.segments.extend(
            starts
                .chunks_exact(2)
                .zip(lens)
                .map(|(pair, &len)| Segment {
                    query: position(pair[0], query_offset),
                    subject: position(pair[1], subject_offset),
                    len: len as u32,
                }),
        );
        self.hsps.push((subject, self.segments.len()));
    }

    /// Each HSP's subject and segments, in order.
    fn iter(&self) -> impl Iterator<Item = (u32, &[Segment])> {
        let mut start = 0;
        self.hsps.iter().map(move |&(subject, end)| {
            let segments = &self.segments[start..end];
            start = end;
            (subject, segments)
        })
    }
}

/// The merged Dense-seg of a section (`x_CreateDenseg`): the sequence of each row (the query
/// first), each segment's length and each segment's starts (`starts[segment * rows + row]`, `-1`
/// where the row has no letters).
pub(crate) struct Denseg {
    ids: Vec<SeqKey>,
    lens: Vec<u32>,
    starts: Vec<i32>,
}

/// What a program gives the section about its sequences.
pub(crate) trait Sequences {
    /// Appends the letters `from..=to` (0-based) of `seq` as the report shows them
    /// (`CAlnVec::GetSeqString`: the Bioseq's IUPAC letters, upper case).
    fn letters(&self, seq: SeqKey, from: u32, to: u32, out: &mut Vec<u8>);

    /// Appends the label of `seq`'s rows (`x_FillSeqid` without `eShowBlastStyleId`).
    fn label(&self, seq: SeqKey, out: &mut Vec<u8>);

    /// Appends the masks of `seq`: the query's `TMaskedQueryRegions`, a subject's of
    /// `GetSubjectMasks` (inclusive, 0-based, in the records; all of frame 0).
    fn masks(&self, seq: SeqKey, out: &mut Vec<(u32, u32)>);
}

/// What a program's report gives the sections of its queries: [`Sequences`] by the query's and
/// the subjects' indices in the report (a subject's masks are those of its HSP list with the
/// query, `GetSubjectMasks` of the query's results).
pub(crate) trait ReportSequences: Sync {
    fn query_letters(&self, q_idx: usize, from: u32, to: u32, out: &mut Vec<u8>);
    fn subject_letters(&self, s_idx: u32, from: u32, to: u32, out: &mut Vec<u8>);
    fn query_label(&self, q_idx: usize, out: &mut Vec<u8>);
    fn subject_label(&self, s_idx: u32, out: &mut Vec<u8>);
    fn query_masks(&self, q_idx: usize, out: &mut Vec<(u32, u32)>);
    fn subject_masks(&self, q_idx: usize, s_idx: u32, out: &mut Vec<(u32, u32)>);
}

/// The [`Sequences`] of query `q_idx`'s section, a subject's key being its index.
pub(crate) struct OfQuery<'a> {
    pub(crate) sequences: &'a dyn ReportSequences,
    pub(crate) q_idx: usize,
}

impl Sequences for OfQuery<'_> {
    fn letters(&self, seq: SeqKey, from: u32, to: u32, out: &mut Vec<u8>) {
        match seq {
            SeqKey::Query => self.sequences.query_letters(self.q_idx, from, to, out),
            SeqKey::Subject(s_idx) => self.sequences.subject_letters(s_idx, from, to, out),
        }
    }

    fn label(&self, seq: SeqKey, out: &mut Vec<u8>) {
        match seq {
            SeqKey::Query => self.sequences.query_label(self.q_idx, out),
            SeqKey::Subject(s_idx) => self.sequences.subject_label(s_idx, out),
        }
    }

    fn masks(&self, seq: SeqKey, out: &mut Vec<(u32, u32)>) {
        match seq {
            SeqKey::Query => self.sequences.query_masks(self.q_idx, out),
            SeqKey::Subject(s_idx) => self.sequences.subject_masks(self.q_idx, s_idx, out),
        }
    }
}

/// The label of the rows of a Bioseq whose Seq-ids are `ids` (`x_FillSeqid` without
/// `eShowBlastStyleId`): with `-show_gis` (`eShowGi`) its GI when it has one, else the content of
/// its `WorstRank` best id with its version. A Bioseq of a local id `Query_N` or `Subject_N` (no
/// `-parse_deflines`) shows that string, one of a database without parsed Seq-ids (`gnl|BL_ORD_ID|n`)
/// its ordinal id `n` (the tag of a `general` id): the callers write those.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4545-4572
/// ```c++
///     } else {
///         if(m_AlignOption&eShowGi){
///             TGi gi = ZERO_GI;
///             if(m_AV->GetSeqId(row).Which() == CSeq_id::e_Gi){
///                 gi = m_AV->GetSeqId(row).GetGi();
///             }
///             if(!(gi > ZERO_GI)){
///                 gi = CAlignFormatUtil::GetGiForSeqIdList(m_AV->GetBioseqHandle(row).\
///                                          GetBioseqCore()->GetId());
///             }
///             if(gi > ZERO_GI){
///                 id=NStr::NumericToString(gi);
///             } else {
///                 const CRef<CSeq_id> wid
///                     = FindBestChoice(m_AV->GetBioseqHandle(row).\
///                                      GetBioseqCore()->GetId(),
///                                      CSeq_id::WorstRank);
///                 id = CAlignFormatUtil::GetLabel(wid,CSeq_id::eContent).c_str();
///             }
///         } else {
///             const CRef<CSeq_id> wid
///                 = FindBestChoice(m_AV->GetBioseqHandle(row).\
///                                  GetBioseqCore()->GetId(),
///                                  CSeq_id::WorstRank);
///             id = CAlignFormatUtil::GetLabel(wid,CSeq_id::eContent).c_str();
///         }
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:919-930
/// ```c++
/// string CAlignFormatUtil::GetLabel(CConstRef<CSeq_id> id,bool with_version)
/// {
///     string retval = "";
///     if (id->Which() == CSeq_id::e_General){
///         const CDbtag& dtg = id->GetGeneral();
///         retval = dtg.GetTag().GetStr();
///     }
///     if(retval == "")
///       retval = id->GetSeqIdString(with_version);
///
///     return retval;
/// }
/// ```
/// (`CSeq_id::eContent` is 1, passed as `with_version` true; `seqid::label::format_label`.)
pub(crate) fn fill_seqid(ids: &[crate::seqid::ParsedSeqId], show_gis: bool, out: &mut Vec<u8>) {
    if show_gis {
        let gi = crate::seqid::label::first_gi(ids);
        if gi > 0 {
            out.extend_from_slice(gi.to_string().as_bytes());
            return;
        }
    }
    if let Some(best) = crate::seqid::rank::find_best(ids, crate::seqid::rank::text_score) {
        crate::seqid::label::format_label(best, true, out);
    }
}

/// Formats 1 to 4: `eShowIdentity` (`.` for a letter equal to the query's: 1 and 3) and
/// `eMasterAnchored` (only the query's columns, with the rows' insertions below them: 1 and 2).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct View {
    pub(crate) identity: bool,
    pub(crate) anchored: bool,
}

/// An error that NCBI's code reaches only with input that a BLAST search does not make.
fn internal(what: &str) -> io::Error {
    io::Error::other(format!(
        "LOSAT internal error in the query-anchored report: {what}"
    ))
}

/// Writes the section of one query whose HSPs are `hsps` (none: nothing, as `DisplaySeqalign`
/// returns before its first byte for an empty Seq-align set; `-num_alignments 0`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1911-1918
/// ```c++
/// void CDisplaySeqalign::DisplaySeqalign(CNcbiOstream& out)
/// {
///     CSeq_align_set actual_aln_list;
///     CAlignFormatUtil::ExtractSeqalignSetFromDiscSegs(actual_aln_list,
///                                                      *m_SeqalignSetRef);
///     if (actual_aln_list.Get().empty()){
///         return;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2055-2111
/// ```c++
///         for(int i = 0; i < (int)alnVector.size(); i ++){
///             bool hasAln = false;
///             for(CTypeConstIterator<CSeq_align>
///                     alnRef = ConstBegin(*alnVector[i]); alnRef; ++alnRef){
///                 CTypeConstIterator<CDense_seg> ds = ConstBegin(*alnRef);
///                 //*out2 << *ds;
///                 try{
///                     if (m_AlignOption & eTranslateNucToNucAlignment) {
///                         mix[i]->Add(*ds, CAlnMix::fForceTranslation);
///                     } else {
///                         if (ds->IsSetWidths() &&
///                             ds->GetWidths()[0] == 3 &&
///                             ds->IsSetStrands() &&
///                             ds->GetStrands().front()==eNa_strand_minus){
///                             mix[i]->Add(*ds, CAlnMix::fNegativeStrand);
///                         } else {
///                             mix[i]->Add(*ds, CAlnMix::fPreserveRows);
///                         }
///                     }
///                 } catch (const CException& e){
///                     _TRACE("Warning: " << e.what());
///                     continue;
///                 }
///                 hasAln = true;
///             }
///             if(hasAln){
///                 //    *out2<<*alnVector[i];
///                 mix[i]->Merge(CAlnMix::fMinGap
///                               | CAlnMix::fQuerySeqMergeOnly
///                               | CAlnMix::fFillUnalignedRegions);
///                 //*out2<<mix[i]->GetDenseg();
///             }
///         }
///
///         int numDistinctFrames = 0;
///         for(int i = 0; i < (int)alnVector.size(); i ++){
///             if(!alnVector[i]->Get().empty()){
///                 numDistinctFrames ++;
///             }
///         }
///         out<<"\n";
///         for(int i = 0; i < k_NumFrame; i ++){
///             try{
///                 CRef<CAlnVec> avRef (new CAlnVec (mix[i]->GetDenseg(),
///                                                   m_Scope));
///                 avRef->SetAaCoding(CSeq_data::e_Ncbieaa);
///                 avRef->SetGenCode(m_SlaveGeneticCode);
///                 avRef->SetGenCode(m_MasterGeneticCode, 0);
///                 m_AV = avRef;
///
///                 if(numDistinctFrames > 1){
///                     out << "For reading frame " << k_FrameConversion[i]
///                         << " of query sequence:\n\n";
///                 }
///                 x_DisplayAlnvec(out);
///             } catch (CException e){
///                 continue;
///             }
///         }
/// ```
/// Every HSP here is a Dense-seg of the first group, added with `fPreserveRows` (which throws
/// nothing for it); the other groups are empty, so `GetDenseg` throws for them and they are
/// skipped, and no frame line is written.
pub(crate) fn write_section<W: Write + ?Sized>(
    writer: &mut W,
    view: View,
    line_length: usize,
    hsps: &Hsps,
    sequences: &dyn Sequences,
) -> io::Result<()> {
    if hsps.is_empty() {
        return Ok(());
    }
    writer.write_all(b"\n")?;
    let denseg = aln_mix::merge(hsps).map_err(|error| internal(error.0))?;
    let mut map = AlnMap::new(&denseg);
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1206-1208
    // ```c++
    //     if(m_AlignOption & eMasterAnchored){
    //         m_AV->SetAnchor(0);
    //     }
    // ```
    if view.anchored {
        map.set_anchor(0);
    }
    let rows = prepare_rows(&map, line_length, sequences);
    display_rows(writer, &map, view, line_length, rows, sequences)
}

/// `SAlnRowInfo`: what `x_PrepareRowData` makes of each row.
struct RowData {
    sequence: Vec<Vec<u8>>,
    strings: Vec<aln_vec::RowString>,
    labels: Vec<Vec<u8>>,
    row_rng: Vec<(i64, i64)>,
    /// `masked_regions`: each row's masks as alignment ranges, by the row's sequence (`masks`).
    masked: Vec<usize>,
    masks: Vec<Vec<(i64, i64)>>,
    max_id_len: usize,
    max_start_len: usize,
}

/// `x_PrepareRowData` for these formats: no features, no statistics, no sequence labels.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1242-1250
/// ```c++
///     // Mask locations for queries (first elem) and subjects (all other rows)
///     vector<TSAlnSeqlocInfoList> masked_regions(rowNum);
///     x_FillLocList(masked_regions[0], m_Seqloc);
///
///     for (int row = 1; row < rowNum; row++) {
///         const CSeq_id& id = m_AV->GetSeqId(row);
///         x_FillLocList(masked_regions[row], &m_SubjectMasks[id]);
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1269-1280
/// ```c++
///         rowRng[row] = m_AV->GetSeqAlnRange(row);
///         frame[row] = (m_AV->GetWidth(row) == 3 ?
///                       ...
///                                  m_AV->GetSeqId(row), m_Scope) : 0);
///         //make sequence
///         m_AV->GetWholeAlnSeqString(row, sequence[row], &insertAlnStart[row],
///                                    &insertStart[row], &insertLength[row],
///                                    (int)m_LineLen, &seqStarts[row], &seqStops[row]);
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1366-1370
/// ```c++
///         //make id
///         x_FillSeqid(seqidArray[row], row);
///         maxIdLen=max<size_t>(seqidArray[row].size(), maxIdLen);
///         size_t maxCood=max<size_t>(m_AV->GetSeqStart(row), m_AV->GetSeqStop(row));
///         maxStartLen = max<size_t>(NStr::SizetToString(maxCood).size(), maxStartLen);
/// ```
/// (`maxStartLen` is the width of the largest 0-based position, which the 1-based numbers
/// printed can exceed by one digit.)
fn prepare_rows(map: &AlnMap<'_>, line_length: usize, sequences: &dyn Sequences) -> RowData {
    let num_rows = map.num_rows();
    let mut letters = |seq: SeqKey, from: u32, to: u32, out: &mut Vec<u8>| {
        sequences.letters(seq, from, to, out);
    };
    let mut data = RowData {
        sequence: Vec::with_capacity(num_rows),
        strings: Vec::with_capacity(num_rows),
        labels: Vec::with_capacity(num_rows),
        row_rng: Vec::with_capacity(num_rows),
        masked: Vec::with_capacity(num_rows),
        masks: Vec::new(),
        max_id_len: 0,
        max_start_len: 0,
    };
    let mut keys: Vec<SeqKey> = Vec::new();
    let mut intervals = Vec::new();
    for row in 0..num_rows {
        let key = map.seq_key(row);
        data.row_rng.push(map.seq_aln_range(row));
        let mut string = map.whole_aln_seq_string(row, line_length as u32, &mut letters);
        data.sequence.push(std::mem::take(&mut string.buffer));
        data.strings.push(string);
        let mut label = Vec::new();
        sequences.label(key, &mut label);
        data.max_id_len = data.max_id_len.max(label.len());
        data.labels.push(label);
        let max_cood = map.seq_start(row).max(map.seq_stop(row));
        data.max_start_len = data.max_start_len.max(max_cood.to_string().len());
        // The rows of one sequence have one list (`x_FillLocList` of its masks).
        let index = match keys.iter().position(|&known| known == key) {
            Some(index) => index,
            None => {
                intervals.clear();
                sequences.masks(key, &mut intervals);
                keys.push(key);
                data.masks.push(fill_loc_list(map, key, &intervals));
                keys.len() - 1
            }
        };
        data.masked.push(index);
    }
    data
}

/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4291-4338
/// ```c++
/// void
/// CDisplaySeqalign::x_FillLocList(TSAlnSeqlocInfoList& loc_list,
///                         const list< CRef<CSeqLocInfo> >* masks) const
/// {
///     if ( !masks ) {
///         return;
///     }
///
///     ITERATE(TMaskedQueryRegions, iter, *masks) {
///         CRef<SAlnSeqlocInfo> alnloc(new SAlnSeqlocInfo);
///         bool has_valid_loc = false;
///         for (int i=0; i<m_AV->GetNumRows(); i++){
///         	const CSeq_interval& interval = (*iter)->GetInterval();
///             TSeqRange loc_range(interval.GetFrom(), interval.GetTo());
///             if(interval.GetId().Match(m_AV->GetSeqId(i)) &&
///                m_AV->GetSeqRange(i).IntersectingWith(loc_range)){
///                 int actualAlnStart = 0, actualAlnStop = 0;
///                 if(m_AV->IsPositiveStrand(i)){
///                     actualAlnStart =
///                         m_AV->GetAlnPosFromSeqPos(i,
///                                                   interval.GetFrom(),
///                                                           CAlnMap::eBackwards, true);
///                     actualAlnStop =
///                         m_AV->GetAlnPosFromSeqPos(i,
///                                                   interval.GetTo(),
///                                                   CAlnMap::eBackwards, true);
///                 } else {
///                     ...
///                 }
///                 alnloc->aln_range.Set(actualAlnStart, actualAlnStop);
///                 has_valid_loc = true;
///                 break;
///             }
///         }
///         if (has_valid_loc) {
///             alnloc->seqloc = *iter;
///             loc_list.push_back(alnloc);
///         }
///     }
/// }
/// ```
/// A mask goes to the columns of the first row of its sequence whose letters it meets, and
/// `x_OutputSeq` puts those columns in lower case in every row of that sequence: a subject's mask
/// spreads to the same columns of its other HSPs' rows (PD-LOSAT-NCBI-DEFECTS rule 2, DECISIONS
/// T3-6).
fn fill_loc_list(map: &AlnMap<'_>, key: SeqKey, masks: &[(u32, u32)]) -> Vec<(i64, i64)> {
    let mut list = Vec::new();
    for &(from, to) in masks {
        for row in 0..map.num_rows() {
            if map.seq_key(row) == key && map.seq_start(row) <= to && from <= map.seq_stop(row) {
                list.push((
                    i64::from(map.aln_pos_from_seq_pos(row, from)),
                    i64::from(map.aln_pos_from_seq_pos(row, to)),
                ));
                break;
            }
        }
    }
    list
}

/// `x_DisplayRowData`: the lines of `line_length` columns.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1435-1451
/// ```c++
/// void CDisplaySeqalign::x_DisplayRowData(SAlnRowInfo *alnRoInfo,CNcbiOstream& out)
/// {
///     size_t aln_stop=m_AV->GetAlnStop();
///     int rowNum = alnRoInfo->rowNum;
///     vector<int> prev_stop(rowNum);
///     ...
///     for(int j=0; j<=(int)aln_stop; j+=(int)m_LineLen){
///         ...
///             string rowdata = x_DisplayRowDataSet(alnRoInfo,j, prev_stop);
///             ...
///             out << rowdata;
/// ```
/// (`j` stays below 2^31: an alignment of more than 2^30 columns would make NCBI's `int` wrap.)
fn display_rows<W: Write + ?Sized>(
    writer: &mut W,
    map: &AlnMap<'_>,
    view: View,
    line_length: usize,
    mut rows: RowData,
    sequences: &dyn Sequences,
) -> io::Result<()> {
    let aln_stop = i64::from(map.aln_stop());
    let mut prev_stop = vec![0i32; map.num_rows()];
    let mut out = Vec::new();
    let mut j = 0i64;
    while j <= aln_stop {
        out.clear();
        display_row_data_set(
            &mut out,
            map,
            view,
            line_length,
            &mut rows,
            j,
            &mut prev_stop,
            sequences,
        )?;
        writer.write_all(&out)?;
        j += line_length as i64;
    }
    Ok(())
}

/// Appends `n` spaces (`CAlignFormatUtil::AddSpace`).
fn add_space(out: &mut Vec<u8>, n: usize) {
    out.extend(std::iter::repeat_n(b' ', n));
}

/// `x_DisplayRowDataSet`: one line of every row that has letters on it, then a blank line.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1477-1554
/// ```c++
///     //output according to aln coordinates
///     if(aln_stop-aln_start+1<m_LineLen) {
///         actualLineLen=aln_stop-aln_start+1;
///     } else {
///         actualLineLen=m_LineLen;
///     }
///     CAlnMap::TSignedRange curRange(aln_start, aln_start+(int)actualLineLen-1);
///     alnRoInfo->currPrintSegment = aln_start;
///     alnRoInfo->currActualLineLen = actualLineLen;
///     alnRoInfo->currRange = curRange;
///     //here is each row
///     for (int row=0; row<rowNum; row++) {
///         bool hasSequence = true;
///         if (!(m_AlignOption & eShowGapOnlyLines)) {
///             hasSequence = curRange.IntersectingWith(alnRoInfo->rowRng[row]);
///         }
///         //only output rows that have sequence
///         if (hasSequence){
///             int end = alnRoInfo->seqStops[row].front() + 1;
///             bool has_mismatch = false;
///             //change the alignment line to identity style
///             if (row>0 && m_AlignOption & eShowIdentity){//check usage - pairwise - only
///                 for (int index = aln_start; index < aln_start + (int)actualLineLen &&
///                          index < (int)alnRoInfo->sequence[row].size(); index ++){
///                     if (alnRoInfo->sequence[row][index] == alnRoInfo->sequence[0][index] &&
///                         isalpha((unsigned char) alnRoInfo->sequence[row][index])) {
///                         alnRoInfo->sequence[row][index] = k_IdentityChar;
///                     } else if (!has_mismatch) {
///                         has_mismatch = true;
///                     }
///                 }
///             }
///             //feature for query
///             if(row == 0){
///                 x_PrintFeatures(alnRoInfo, row, master_feat_str, out);
///             }
///             if((m_AlignOption & eMergeAlign) || (m_AlignOption & eHyperLinkMasterSeqid) || (m_AlignOption & eHyperLinkSlaveSeqid)) {
///                 x_DisplaySequenceIDForQueryAnchored(alnRoInfo,row,out);
///             }
///             else {
///                 ...
///             }
///             //print out sequence line
///             x_DisplaySequenceLine(alnRoInfo, row, prev_stop[row], out);
///             if(m_AlignOption & eMasterAnchored){
///                 //inserts for anchored view
///                 x_DisplayInsertsForQueryAnchored(alnRoInfo,row,out);
///             }
///             ...
///             prev_stop[row] = end;
///         }
///         if(!alnRoInfo->seqStarts[row].empty()){ //shouldn't need this check
///                 alnRoInfo->seqStarts[row].pop_front();
///         }
///         if(!alnRoInfo->seqStops[row].empty()){
///                 alnRoInfo->seqStops[row].pop_front();
///         }
///     }//end of displaying rows
///     out<<"\n";
/// ```
/// The identity test reads the letters before any mask puts them in lower case (`x_OutputSeq`
/// masks a copy), so a masked identical letter is a `.`. The label is the whole id
/// (`x_DisplaySequenceIDForQueryAnchored` without HTML, alignment statistics or sequence
/// property labels).
#[allow(clippy::too_many_arguments)]
fn display_row_data_set(
    out: &mut Vec<u8>,
    map: &AlnMap<'_>,
    view: View,
    line_length: usize,
    rows: &mut RowData,
    aln_start: i64,
    prev_stop: &mut [i32],
    sequences: &dyn Sequences,
) -> io::Result<()> {
    let aln_stop = i64::from(map.aln_stop());
    let actual_line_len = if aln_stop - aln_start + 1 < line_length as i64 {
        (aln_stop - aln_start + 1) as usize
    } else {
        line_length
    };
    let cur_range = (aln_start, aln_start + actual_line_len as i64 - 1);
    for row in 0..map.num_rows() {
        let (from, to) = rows.row_rng[row];
        let has_sequence = cur_range.0 <= to && from <= cur_range.1;
        if has_sequence {
            let front = |list: &VecDeque<u32>| -> io::Result<i32> {
                list.front()
                    .map(|&value| value.wrapping_add(1) as i32)
                    .ok_or_else(|| internal("a row has no coordinates for a line"))
            };
            let end = front(&rows.strings[row].rights)?;
            let start = front(&rows.strings[row].lefts)?;
            if row > 0 && view.identity {
                let (query, others) = rows.sequence.split_at_mut(1);
                let letters = &mut others[row - 1];
                let stop = (aln_start as usize + actual_line_len).min(letters.len());
                for index in aln_start as usize..stop {
                    if letters[index] == query[0][index] && letters[index].is_ascii_alphabetic() {
                        letters[index] = b'.';
                    }
                }
            }
            out.extend_from_slice(&rows.labels[row]);
            display_sequence_line(
                out,
                rows,
                row,
                start,
                end,
                aln_start,
                actual_line_len,
                prev_stop[row],
            );
            if view.anchored {
                display_inserts(
                    out,
                    map,
                    rows,
                    row,
                    aln_start,
                    actual_line_len,
                    line_length,
                    sequences,
                )?;
            }
            prev_stop[row] = end;
        }
        rows.strings[row].lefts.pop_front();
        rows.strings[row].rights.pop_front();
    }
    out.push(b'\n');
    Ok(())
}

/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1594-1627
/// ```c++
/// void CDisplaySeqalign::x_DisplaySequenceLine(SAlnRowInfo *alnRoInfo, int row, int prev_stop, CNcbiOstrstream &out)
/// {
///     size_t startLen=0;
///     int start = alnRoInfo->seqStarts[row].front() + 1;  //+1 for 1 based
///     int end = alnRoInfo->seqStops[row].front() + 1;
///     int j = alnRoInfo->currPrintSegment;
///     int actualLineLen = alnRoInfo->currActualLineLen;
///     //print out sequence line
///     //adjust space between id and start
///     CAlignFormatUtil::AddSpace(out, alnRoInfo->maxIdLen-alnRoInfo->seqidArray[row].size() + k_IdStartMargin);
///     //not to display start and stop number for empty row
///     if ((j > 0 && end == prev_stop)
///         || (j == 0 && start == 1 && end == 1)) {
///         startLen = 0;
///     } else {
///         out << start;
///         startLen=NStr::IntToString(start).size();
///     }
///
///     CAlignFormatUtil::AddSpace(out, alnRoInfo->maxStartLen-startLen + k_StartSequenceMargin);
///     x_OutputSeq(alnRoInfo->sequence[row], m_AV->GetSeqId(row), j,
///                 (int)actualLineLen, alnRoInfo->frame[row], row,
///                 (row > 0 && alnRoInfo->colorMismatch)?true:false,
///                 alnRoInfo->masked_regions[row], out);
///     CAlignFormatUtil::AddSpace(out, k_SeqStopMargin);
///
///      //not to display stop number for empty row in the middle
///     if (!(j > 0 && end == prev_stop)
///         && !(j == 0 && start == 1 && end == 1)) {
///         out << end;
///     }
///     out<<"\n";
/// }
/// ```
/// (`k_IdStartMargin`, `k_StartSequenceMargin` and `k_SeqStopMargin` are 2, showalign.cpp:96-98;
/// the `size_t` sum `maxStartLen - startLen + 2` is at least 1.)
#[allow(clippy::too_many_arguments)]
fn display_sequence_line(
    out: &mut Vec<u8>,
    rows: &RowData,
    row: usize,
    start: i32,
    end: i32,
    j: i64,
    actual_line_len: usize,
    prev_stop: i32,
) {
    add_space(out, rows.max_id_len - rows.labels[row].len() + 2);
    let empty = (j > 0 && end == prev_stop) || (j == 0 && start == 1 && end == 1);
    let start_len = if empty {
        0
    } else {
        let text = start.to_string();
        out.extend_from_slice(text.as_bytes());
        text.len()
    };
    add_space(out, (rows.max_start_len + 2).wrapping_sub(start_len));
    output_seq(
        out,
        &rows.sequence[row],
        j as usize,
        actual_line_len,
        &rows.masks[rows.masked[row]],
    );
    add_space(out, 2);
    if !empty {
        out.extend_from_slice(end.to_string().as_bytes());
    }
    out.push(b'\n');
}

/// `x_OutputSeq` of a row line (no HTML): the columns `start..start + len` of `sequence`, those
/// of the row's masks in lower case (`m_SeqLocChar` is `eLowerCase`, blast_format.cpp:1581).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2486-2523
/// ```c++
///     _ASSERT((int)sequence.size() > start);
///     list<CRange<int> > actualSeqloc;
///     string actualSeq = sequence.substr(start, len);
///
///     if(id.Which() != CSeq_id::e_not_set){
///         /*only do this for sequence but not for others like middle line,
///           features*/
///         ITERATE(TSAlnSeqlocInfoList, iter, loc_list) {
///             int from=(*iter)->aln_range.GetFrom();
///             int to=(*iter)->aln_range.GetTo();
///             int locFrame = (*iter)->seqloc->GetFrame();
///             if(id.Match((*iter)->seqloc->GetInterval().GetId())
///                && locFrame == frame){
///                 bool isFirstChar = true;
///                 CRange<int> eachSeqloc(0, 0);
///                 //go through each residule and mask it
///                 for (int i=max<int>(from, start);
///                      i<=min<int>(to, start+len -1); i++){
///                     ...
///                     } else if (m_SeqLocChar==eLowerCase){
///                         actualSeq[i-start]=tolower((unsigned char) actualSeq[i-start]);
///                     }
/// ```
/// (The row's masks are those of its own sequence, all of frame 0, as is the row.)
fn output_seq(out: &mut Vec<u8>, sequence: &[u8], start: usize, len: usize, masks: &[(i64, i64)]) {
    let end = (start + len).min(sequence.len());
    let first = out.len();
    out.extend_from_slice(&sequence[start.min(end)..end]);
    for &(from, to) in masks {
        let lo = from.max(start as i64);
        let hi = to.min(start as i64 + len as i64 - 1);
        let mut column = lo;
        while column <= hi {
            if let Some(letter) = out.get_mut(first + (column - start as i64) as usize) {
                letter.make_ascii_lowercase();
            }
            column += 1;
        }
    }
}

/// One insertion of a row on a line (`SInsertInformation`): the column after which it goes
/// (`aln_start`, the anchor's columns before it less one), its first letter and its length.
#[derive(Clone, Copy)]
struct Insert {
    aln_start: i64,
    seq_start: u32,
    insert_len: u32,
}

/// `x_DisplayInsertsForQueryAnchored`: a line with `\` after each insertion point of the row on
/// this line, then for each round of `x_DoFills` a line of bars and a line of letters.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1629-1677
/// ```c++
///     list<string> inserts;
///     string insertPosString;  //the one with "\" to indicate insert
///     TSInsertInformationList insertList;
///     int j = alnRoInfo->currPrintSegment;
///     CAlnMap::TSignedRange curRange = alnRoInfo->currRange;
///     x_GetInserts(insertList, alnRoInfo->insertAlnStart[row],
///                                  alnRoInfo->insertStart[row], alnRoInfo->insertLength[row],
///                                  j + (int)m_LineLen);
///     x_FillInserts(row, curRange, j, inserts, insertPosString, insertList);
///     bool insertAlready = false;
///     for(list<string>::iterator iter = inserts.begin();
///         iter != inserts.end(); iter ++){
///         if(!insertAlready){
///             ...
///             int base_margin = alnRoInfo->maxIdLen + k_IdStartMargin + alnRoInfo->maxStartLen + k_StartSequenceMargin;
///             ...
///             CAlignFormatUtil::AddSpace(out, base_margin);
///             out << insertPosString<<"\n";
///         }
///         ...
///         int base_margin = alnRoInfo->maxIdLen + k_IdStartMargin + alnRoInfo->maxStartLen + k_StartSequenceMargin;
///         ...
///         CAlignFormatUtil::AddSpace(out, base_margin);
///         out<<*iter<<"\n";
///         insertAlready = true;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3099-3120
/// ```c++
///     while(!insert_aln_start.empty()
///           && (int)insert_aln_start.front() < line_aln_stop){
///         CRef<SInsertInformation> insert(new SInsertInformation);
///         insert->aln_start = insert_aln_start.front() - 1;
///         //Need to minus one as we are inserting after this position
///         insert->seq_start = insert_seq_start.front();
///         insert->insert_len = insert_length.front();
///         insert_list.push_back(insert);
///         insert_aln_start.pop_front();
///         insert_seq_start.pop_front();
///         insert_length.pop_front();
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3076-3096
/// ```c++
///     string line(aln_range.GetLength(), ' ');
///
///     ITERATE(TSInsertInformationList, iter, insert_list){
///         int from = (*iter)->aln_start;
///         line[from - aln_start + 1] = '\\';
///     }
///     insert_pos_string = line;
///     //this is the line with "\" right after each insert position
///
///     //here fills the insert sequence
///     x_DoFills(row, aln_range, aln_start, insert_list, inserts);
/// ```
#[allow(clippy::too_many_arguments)]
fn display_inserts(
    out: &mut Vec<u8>,
    map: &AlnMap<'_>,
    rows: &mut RowData,
    row: usize,
    aln_start: i64,
    actual_line_len: usize,
    line_length: usize,
    sequences: &dyn Sequences,
) -> io::Result<()> {
    let line_aln_stop = aln_start + line_length as i64;
    let mut insert_list = Vec::new();
    let queue = &mut rows.strings[row].inserts;
    while let Some(&(aln, seq_start, insert_len)) = queue.front() {
        if i64::from(aln as i32) >= line_aln_stop {
            break;
        }
        insert_list.push(Insert {
            aln_start: i64::from(aln.wrapping_sub(1) as i32),
            seq_start,
            insert_len,
        });
        queue.pop_front();
    }
    let mut line = vec![b' '; actual_line_len];
    for insert in &insert_list {
        *column(&mut line, insert.aln_start - aln_start + 1)? = b'\\';
    }
    let mut inserts = Vec::new();
    do_fills(
        map.seq_key(row),
        actual_line_len,
        aln_start,
        insert_list,
        &mut inserts,
        sequences,
    )?;
    let base_margin = rows.max_id_len + 2 + rows.max_start_len + 2;
    for (index, insert) in inserts.iter().enumerate() {
        if index == 0 {
            add_space(out, base_margin);
            out.extend_from_slice(&line);
            out.push(b'\n');
        }
        add_space(out, base_margin);
        out.extend_from_slice(insert);
        out.push(b'\n');
    }
    Ok(())
}

/// The byte of `line` at `index`: NCBI writes the marks of an insertion inside the line (an
/// insertion comes before a query column of the line).
fn column(line: &mut [u8], index: i64) -> io::Result<&mut u8> {
    usize::try_from(index)
        .ok()
        .and_then(|index| line.get_mut(index))
        .ok_or_else(|| internal("an insertion outside its line"))
}

/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3031-3073
/// ```c++
/// void CDisplaySeqalign::x_DoFills(int row, CAlnMap::TSignedRange& aln_range,
///                                  int  aln_start,
///                                  TSInsertInformationList& insert_list,
///                                  list<string>& inserts) const {
///     if(!insert_list.empty()){
///         string bar(aln_range.GetLength(), ' ');
///
///         string seq;
///         TSInsertInformationList leftOverInsertList;
///         bool isFirstInsert = true;
///         int curInsertAlnStart = 0;
///         int prvsInsertAlnEnd = 0;
///
///         //go through each insert and fills the seq if it can
///         //be filled on the same line.  If not, go to the next line
///         NON_CONST_ITERATE(TSInsertInformationList, iter, insert_list) {
///             curInsertAlnStart = (*iter)->aln_start;
///             //always fill the first insert.  Also fill if there is enough space
///             if(isFirstInsert || curInsertAlnStart - prvsInsertAlnEnd >= 1){
///                 bar[curInsertAlnStart-aln_start+1] = '|';
///                 int seqStart = (*iter)->seq_start;
///                 int seqEnd = seqStart + (*iter)->insert_len - 1;
///                 string newInsert;
///                 newInsert = m_AV->GetSeqString(newInsert, row, seqStart,
///                                                seqEnd);
///                 prvsInsertAlnEnd = s_AdjustInsert(seq, newInsert,
///                                                   curInsertAlnStart, aln_start);
///                 isFirstInsert = false;
///             } else { //if no space, save the chunk and go to next line
///                 bar[curInsertAlnStart-aln_start+1] = '|';
///                 //indicate insert goes to the next line
///                 prvsInsertAlnEnd += x_AddBar(seq, curInsertAlnStart, aln_start);
///                 //May need to add a bar after the current insert sequence
///                 //to indicate insert goes to the next line.
///                 leftOverInsertList.push_back(*iter);
///             }
///         }
///         //save current insert.  Note that each insert has a bar and sequence
///         //below it
///         inserts.push_back(bar);
///         inserts.push_back(seq);
///         //here recursively fill the chunk that don't have enough space
///         x_DoFills(row, aln_range, aln_start, leftOverInsertList, inserts);
///     }
///
/// }
/// ```
/// (The recursion as a loop.) The letters of an insertion are not masked.
fn do_fills(
    key: SeqKey,
    aln_range_len: usize,
    aln_start: i64,
    mut insert_list: Vec<Insert>,
    inserts: &mut Vec<Vec<u8>>,
    sequences: &dyn Sequences,
) -> io::Result<()> {
    let mut new_insert = Vec::new();
    while !insert_list.is_empty() {
        let mut bar = vec![b' '; aln_range_len];
        let mut seq = Vec::new();
        let mut left_over = Vec::new();
        let mut is_first_insert = true;
        let mut prvs_insert_aln_end = 0i64;
        for insert in insert_list {
            let cur_insert_aln_start = insert.aln_start;
            *column(&mut bar, cur_insert_aln_start - aln_start + 1)? = b'|';
            if is_first_insert || cur_insert_aln_start - prvs_insert_aln_end >= 1 {
                new_insert.clear();
                let seq_end = insert.seq_start + insert.insert_len - 1;
                sequences.letters(key, insert.seq_start, seq_end, &mut new_insert);
                prvs_insert_aln_end =
                    adjust_insert(&mut seq, &new_insert, cur_insert_aln_start, aln_start);
                is_first_insert = false;
            } else {
                prvs_insert_aln_end += add_bar(&mut seq, cur_insert_aln_start, aln_start);
                left_over.push(insert);
            }
        }
        inserts.push(bar);
        inserts.push(seq);
        insert_list = left_over;
    }
    Ok(())
}

/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2980-2992
/// ```c++
/// static int x_AddBar(string& seq, int insert_alnpos, int aln_start){
///     int end = (int)seq.size() -1 ;
///     int barPos = insert_alnpos - aln_start + 1;
///     string addOn;
///     if(barPos - end > 1){
///         string spacer(barPos - end - 1, ' ');
///         addOn += spacer + "|";
///     } else if (barPos - end == 1){
///         addOn += "|";
///     }
///     seq += addOn;
///     return max<int>((barPos - end), 0);
/// }
/// ```
fn add_bar(seq: &mut Vec<u8>, insert_alnpos: i64, aln_start: i64) -> i64 {
    let end = seq.len() as i64 - 1;
    let bar_pos = insert_alnpos - aln_start + 1;
    if bar_pos - end > 1 {
        add_space(seq, (bar_pos - end - 1) as usize);
        seq.push(b'|');
    } else if bar_pos - end == 1 {
        seq.push(b'|');
    }
    (bar_pos - end).max(0)
}

/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3003-3029
/// ```c++
/// static int s_AdjustInsert(string& cur_insert, string& new_insert,
///                           int insert_alnpos, int aln_start)
/// {
///     int insertEnd = 0;
///     int curInsertSize = (int)cur_insert.size();
///     int insertLeftSpace = insert_alnpos - aln_start - curInsertSize + 2;
///     //plus2 because insert is put after the position
///     if(curInsertSize > 0){
///         _ASSERT(insertLeftSpace >= 2);
///     }
///     int newInsertSize = (int)new_insert.size();
///     if(insertLeftSpace - newInsertSize >= 1){
///         //can insert with the end position right below the bar
///         string spacer(insertLeftSpace - newInsertSize, ' ');
///         cur_insert += spacer + new_insert;
///
///     } else { //Need to insert beyond the insert postion
///         if(curInsertSize > 0){
///             cur_insert += " " + new_insert;
///         } else {  //can insert right at the firt position
///             cur_insert += new_insert;
///         }
///     }
///     insertEnd = aln_start + (int)cur_insert.size() -1 ; //-1 back to string position
///     return insertEnd;
/// }
/// ```
fn adjust_insert(
    cur_insert: &mut Vec<u8>,
    new_insert: &[u8],
    insert_alnpos: i64,
    aln_start: i64,
) -> i64 {
    let cur_insert_size = cur_insert.len() as i64;
    let insert_left_space = insert_alnpos - aln_start - cur_insert_size + 2;
    let new_insert_size = new_insert.len() as i64;
    if insert_left_space - new_insert_size >= 1 {
        add_space(cur_insert, (insert_left_space - new_insert_size) as usize);
        cur_insert.extend_from_slice(new_insert);
    } else {
        if cur_insert_size > 0 {
            cur_insert.push(b' ');
        }
        cur_insert.extend_from_slice(new_insert);
    }
    aln_start + cur_insert.len() as i64 - 1
}

#[cfg(test)]
mod tests;
