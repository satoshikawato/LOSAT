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
//!
//! What the programs' Seq-aligns reach (lines RG T3a, T3b and T3c): Dense-segs of two rows
//! (BLASTN's minus-strand HSPs with the query on the minus strand; the Dense-diag of an ungapped
//! BLASTN HSP, `CreateDensegFromDendiag`, is the Dense-seg of its one segment), and the Std-segs
//! of the translated programs as `CreateDensegFromStdseg` makes them: TBLASTN's (a protein query
//! row of width 1, a nucleotide subject row of width 3 in its frame), BLASTX's (a nucleotide query
//! row of width 3, a protein subject row) and TBLASTX's (no widths: two nucleotide rows of the same
//! length, which `CAlnMix::Add` with `fForceTranslation` gives widths of 3, `x_ExtendDSWithWidths`).
//! A translated query's HSPs go into six groups by the query's frame (`s_GetStdsegMasterFrame`),
//! each merged and written on its own, under `For reading frame <f> of query sequence:` when more
//! than one group has HSPs (`write_section`). BLASTN's, BLASTP's and TBLASTN's are all in the
//! first group, added with `fPreserveRows`; TBLASTX's with `fForceTranslation`, and BLASTX's of a
//! minus frame with `fNegativeStrand`, which `Add` reads as `fForceTranslation` (`AddFlags`).

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

/// The strands of an HSP's Dense-seg rows and the query's frame (`Hsps::push_dense_seg`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct HspStrands {
    /// The query's row is on the minus strand: a minus-strand BLASTN HSP (`s_CollectSeqAlignData`
    /// puts the query on the minus strand), a BLASTX or TBLASTX HSP in a minus query frame.
    pub(crate) query_minus: bool,
    /// The subject's row is on the minus strand: a TBLASTN or TBLASTX HSP in a minus subject frame.
    pub(crate) subject_minus: bool,
    /// The query's frame in its record (`s_GetStdsegMasterFrame`: 1 to 3, -1 to -3) when the query
    /// is translated (BLASTX, TBLASTX); 0 otherwise.
    pub(crate) query_frame: i8,
    /// The Dense-seg has no widths: `CreateDensegFromStdseg` of the empty rows of a BLASTX HSP of
    /// length 0 (`serial::seqalign::empty_row_segments`) determines none and resets them
    /// (Seq_align.cpp:906-935). `DisplaySeqalign` puts such a Dense-seg in the first group, the
    /// query's frame +1, whatever its frame, and adds it with `fPreserveRows` (`Hsps::groups`,
    /// `Hsps::add_flags`).
    pub(crate) no_widths: bool,
}

impl HspStrands {
    /// `m_StrandsDiffer` of the HSP's matches (alnmatch.cpp:241-251, cited at `aln_mix::Mix::add`):
    /// one row on the minus strand and the other not (a protein row has no strand, read as plus).
    pub(crate) fn differ(self) -> bool {
        self.query_minus != self.subject_minus
    }
}

/// One HSP of a section (`Hsps::iter`): its subject's key, its Dense-seg's segments and its rows'
/// strands.
#[derive(Clone, Copy)]
pub(crate) struct Hsp<'a> {
    pub(crate) subject: u32,
    pub(crate) segments: &'a [Segment],
    pub(crate) strands: HspStrands,
}

/// The HSPs of one query's section in the order of its Seq-align set after the sort and the
/// prune of the report (`PruneSeqalign`): each with its subject's key and its segments.
pub(crate) struct Hsps {
    segments: Vec<Segment>,
    /// Each HSP's subject, the end of its segments in `segments` and `Hsp::strands`.
    hsps: Vec<(u32, usize, HspStrands)>,
    /// The subject that is the query's own Bioseq (one record read as the query and as a subject
    /// with `-parse_deflines`), which `CAlnMix::Add` keeps on its own row (`fPreserveRows`).
    query_bioseq: Option<u32>,
    /// The width of the query's row and of the subjects' rows: 3 for a nucleotide sequence
    /// translated against proteins (BLASTX's query, TBLASTN's subjects, TBLASTX's both), whose
    /// starts are nucleotides and whose segments' lengths are residues; else 1.
    query_width: u32,
    subject_width: u32,
    /// TBLASTX (`eTranslateNucToNucAlignment`, blast_format.cpp:651): the HSPs go into the groups of
    /// their query frames and are added with `fForceTranslation`.
    translate_nuc_to_nuc: bool,
}

impl Default for Hsps {
    /// The HSPs of a section of untranslated rows (BLASTP, BLASTN).
    fn default() -> Self {
        Self {
            segments: Vec::new(),
            hsps: Vec::new(),
            query_bioseq: None,
            query_width: 1,
            subject_width: 1,
            translate_nuc_to_nuc: false,
        }
    }
}

impl Hsps {
    /// The HSPs of a section whose subjects are nucleotides translated against a protein query
    /// (TBLASTN).
    pub(crate) fn of_translated_subjects() -> Self {
        Self {
            subject_width: 3,
            ..Self::default()
        }
    }

    /// The HSPs of a section whose query is a nucleotide sequence translated against protein
    /// subjects (BLASTX): grouped by the query's frame.
    pub(crate) fn of_translated_query() -> Self {
        Self {
            query_width: 3,
            ..Self::default()
        }
    }

    /// The HSPs of a section of nucleotide sequences translated on both sides (TBLASTX): grouped by
    /// the query's frame and added with `fForceTranslation`.
    pub(crate) fn of_translated_nucleotides() -> Self {
        Self {
            query_width: 3,
            subject_width: 3,
            translate_nuc_to_nuc: true,
            ..Self::default()
        }
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.hsps.is_empty()
    }

    /// Adds an HSP of subject `subject`, its Dense-seg's `starts` (two per segment, the query's
    /// then the subject's, `-1` for a gap, `report::serial::seqalign::Segments::starts`) and
    /// `lens`, each start moved by `query_offset` and `subject_offset` (the remaps of the
    /// ranges, `RemapToQueryLoc` and `s_RemapToSubjectLoc`), and its rows' `strands`.
    pub(crate) fn push_dense_seg(
        &mut self,
        subject: u32,
        starts: &[i32],
        lens: &[i32],
        query_offset: u32,
        subject_offset: u32,
        strands: HspStrands,
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
        self.hsps.push((subject, self.segments.len(), strands));
    }

    /// Marks subject `subject` as the query's own Bioseq (`query_bioseq`).
    pub(crate) fn set_query_bioseq(&mut self, subject: u32) {
        self.query_bioseq = Some(subject);
    }

    /// HSP `index`.
    fn hsp(&self, index: usize) -> Hsp<'_> {
        let (subject, end, strands) = self.hsps[index];
        let start = index.checked_sub(1).map_or(0, |before| self.hsps[before].1);
        Hsp {
            subject,
            segments: &self.segments[start..end],
            strands,
        }
    }

    /// The HSPs of each frame group, by index (`k_FrameConversion`: +1, +2, +3, -1, -2, -3): the
    /// group of the query's frame when the query is translated (BLASTX: the Dense-seg's first
    /// width is 3; TBLASTX: `eTranslateNucToNucAlignment`), else the first group.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2006-2041
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:411-450
    /// The query's first Std-seg holds the HSP's first letter on the query's strand: its frame is
    /// that of the HSP in the query's record (`HspStrands::query_frame`, which the program gives).
    /// A BLASTX Dense-seg without widths (`HspStrands::no_widths`, an HSP of length 0) is in the
    /// first group (the `else` above).
    fn groups(&self) -> [Vec<usize>; K_NUM_FRAME] {
        let mut groups: [Vec<usize>; K_NUM_FRAME] = Default::default();
        let by_frame = self.query_width == 3 || self.translate_nuc_to_nuc;
        for (index, &(_, _, strands)) in self.hsps.iter().enumerate() {
            let group = if !by_frame || (strands.no_widths && !self.translate_nuc_to_nuc) {
                Some(0)
            } else {
                match strands.query_frame {
                    1 => Some(0),
                    2 => Some(1),
                    3 => Some(2),
                    -1 => Some(3),
                    -2 => Some(4),
                    -3 => Some(5),
                    _ => None,
                }
            };
            if let Some(group) = group {
                groups[group].push(index);
            }
        }
        groups
    }

    /// The flag that `DisplaySeqalign` adds HSP `strands`' Dense-seg with.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2061-2072
    /// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnmix.hpp:69-103
    /// `fNegativeStrand` is a merge flag: as the add flag of a BLASTX HSP in a minus query frame
    /// (the Dense-seg's first width is 3 and its first strand minus) it is `fForceTranslation`. A
    /// Dense-seg without widths (`HspStrands::no_widths`) is added with `fPreserveRows`.
    fn add_flags(&self, strands: HspStrands) -> AddFlags {
        if self.translate_nuc_to_nuc
            || (self.query_width == 3 && strands.query_minus && !strands.no_widths)
        {
            AddFlags::ForceTranslation
        } else {
            AddFlags::PreserveRows
        }
    }
}

/// `k_NumFrame`, the frame groups of a section (showalign.cpp:85-86).
const K_NUM_FRAME: usize = 6;

/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:86
const K_FRAME_CONVERSION: [&str; K_NUM_FRAME] = ["+1", "+2", "+3", "-1", "-2", "-3"];

/// The add flags that `DisplaySeqalign` gives `CAlnMix::Add` (`Hsps::add_flags`): `fPreserveRows`
/// (0x04), or `fForceTranslation` (0x02, also as the value of the merge flag `fNegativeStrand`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum AddFlags {
    PreserveRows,
    ForceTranslation,
}

/// The merged Dense-seg of a section (`x_CreateDenseg`): the sequence of each row (the query
/// first), each segment's length and each segment's starts (`starts[segment * rows + row]`, `-1`
/// where the row has no letters), each row's strand (`plus`) and, in a section of translated
/// subjects, each row's width (`widths`, empty otherwise: 1).
pub(crate) struct Denseg {
    ids: Vec<SeqKey>,
    lens: Vec<u32>,
    starts: Vec<i32>,
    plus: Vec<bool>,
    widths: Vec<u32>,
}

/// Outfmt 0's `CAlnVec` of one HSP (`CDisplaySeqalign::x_GetAlnVecForSeqalign`, without an
/// anchor: outfmt 0 has no `eMasterAnchored`): the HSP's own Dense-seg of two rows, the query's
/// and the subject's, each row's alignment range (`GetSeqAlnRange`) and its first and last
/// position on each line of `-line_length` columns (`scrn_lefts` and `scrn_rights` of
/// `GetWholeAlnSeqString`, `aln_vec::AlnMap::screen_coords`), which `report::pairwise` prints.
/// Kept from HSP to HSP to reuse its vectors (`report::pairwise` keeps one per thread).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1823-1861
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:979-992
/// (`CreateTranslatedDensegFromNADenseg`)
/// The HSP's segments are those of its Seq-align (`report::serial::seqalign`: the starts in the
/// records, the lengths in residues, `CreateDensegFromStdseg`'s for a Std-seg) with the strands
/// that `s_CollectSeqAlignData` gives them, the Dense-seg's (`CAlnMap::IsPositiveStrand` reads a
/// row's strand of the first segment, alnmap.hpp:599-603, cited at `AlnMap::width`). TBLASTX's
/// Std-seg of two nucleotide rows becomes one of widths 3 with the lengths in residues
/// (`CreateTranslatedDensegFromNADenseg`), which the segments already have. The Dense-diag of an
/// ungapped BLASTN HSP (`CreateDensegFromDendiag`) is the Dense-seg of its one segment.
pub(crate) struct HspAlnVec {
    ds: Denseg,
    map: aln_vec::MapBuffers,
    rows: [HspRow; 2],
    aln_stop: u32,
}

/// A row of an `HspAlnVec`: `scrn_lefts` and `scrn_rights` (`TSeqPosList`, whose fronts are
/// read and popped a line at a time: the entries of line `k` are the `k`-th, `HspAlnVec::line`),
/// as pairs, and `GetSeqAlnRange`.
#[derive(Default)]
struct HspRow {
    screens: Vec<(u32, u32)>,
    range: (i64, i64),
}

impl Default for HspAlnVec {
    fn default() -> Self {
        Self {
            ds: Denseg {
                // The rows' sequences, which `set` keeps.
                ids: vec![SeqKey::Query, SeqKey::Subject(0)],
                lens: Vec::new(),
                starts: Vec::new(),
                plus: Vec::new(),
                widths: Vec::new(),
            },
            map: aln_vec::MapBuffers::default(),
            rows: Default::default(),
            aln_stop: 0,
        }
    }
}

/// An error that NCBI's outfmt 0 code reaches only with input that a BLAST search does not make.
fn pairwise_internal(what: &str) -> io::Error {
    io::Error::other(format!(
        "LOSAT internal error in the pairwise report: {what}"
    ))
}

impl HspAlnVec {
    /// Makes the `CAlnVec` of the HSP whose Seq-align has the segments `starts` (two per
    /// segment, the query's then the subject's, `-1` for a gap), `lens` (in residues) and
    /// `strands` (two per segment, `report::serial::seqalign::Segments`), its rows of widths
    /// `widths` (the query's and the subject's: 3 for a row of nucleotides translated against
    /// proteins, `CreateDensegFromStdseg`; `[1, 1]` for a Dense-seg without widths), and its rows'
    /// per-line coordinates for lines of `line_length` columns (`(int)m_LineLen`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1281-1283
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1272
    pub(crate) fn set(
        &mut self,
        starts: &[i32],
        lens: &[i32],
        strands: &[crate::report::serial::seqalign::Strand],
        widths: [u32; 2],
        line_length: usize,
    ) -> io::Result<()> {
        use crate::report::serial::seqalign::Strand;
        if starts.len() != 2 * lens.len() || strands.len() != starts.len() || lens.is_empty() {
            return Err(pairwise_internal(
                "an HSP's segments are not a Dense-seg of two rows",
            ));
        }
        let ds = &mut self.ds;
        ds.lens.clear();
        ds.lens.extend(lens.iter().map(|&len| len as u32));
        ds.starts.clear();
        ds.starts.extend_from_slice(starts);
        ds.widths.clear();
        if widths != [1, 1] {
            ds.widths.extend(widths);
        }
        let (mut query, mut subject) = (strands[0], strands[1]);
        if query == Strand::Minus && widths[0] != 3 {
            ds.reverse();
            (query, subject) = (query.reversed(), subject.reversed());
        }
        ds.plus.clear();
        ds.plus
            .extend([query != Strand::Minus, subject != Strand::Minus]);
        let map = AlnMap::new_in(&self.ds, &mut self.map).map_err(pairwise_internal)?;
        for (row, data) in self.rows.iter_mut().enumerate() {
            let screens = &mut data.screens;
            screens.clear();
            map.screen_coords(row, line_length as u32, |left, right| {
                screens.push((left, right));
            })?;
            data.range = map.seq_aln_range(row);
        }
        self.aln_stop = map.aln_stop();
        Ok(())
    }

    /// `GetAlnStop()` of the HSP (`x_DisplayRowData`'s `aln_stop`): `u32::MAX` for an HSP of
    /// length 0 (BLASTX's empty rows), whose lines `j <= (int)aln_stop` are none.
    pub(crate) fn aln_stop(&self) -> u32 {
        self.aln_stop
    }

    /// Row `row` on line `line` (0 for the first), of columns `from..=to`: its first and last
    /// position there (`seqStarts[row].front()`, `seqStops[row].front()`, 0-based) when the row
    /// meets the line (`hasSequence`), `None` otherwise. Every line pops each row's fronts
    /// (`pop_front`), so the fronts of line `line` are the row's `line`-th entries.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1497-1504,1544-1549
    /// (A row has an entry for each line: its screens are at most as wide as the lines.)
    #[inline(always)]
    pub(crate) fn line(
        &self,
        row: usize,
        line: usize,
        from: i64,
        to: i64,
    ) -> io::Result<Option<(u32, u32)>> {
        let data = &self.rows[row];
        let (row_from, row_to) = data.range;
        if !(from <= row_to && row_from <= to) {
            return Ok(None);
        }
        match data.screens.get(line) {
            Some(&screen) => Ok(Some(screen)),
            None => Err(pairwise_internal("a row has no coordinates for a line")),
        }
    }
}

impl Denseg {
    /// `CDense_seg::Reverse`: the segments in the other order, each row on the other strand
    /// (here the rows' strands are the caller's, `HspAlnVec::set`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:661-680
    fn reverse(&mut self) {
        let dim = self.ids.len();
        self.lens.reverse();
        let num_segs = self.lens.len();
        for seg in 0..num_segs / 2 {
            for row in 0..dim {
                self.starts
                    .swap(seg * dim + row, (num_segs - 1 - seg) * dim + row);
            }
        }
    }
}

/// What a program gives the section about its sequences.
pub(crate) trait Sequences {
    /// Appends the letters `from..=to` (0-based, in the records) of `seq`'s plus strand as the
    /// report reads them (the IUPAC letters of the Bioseq's sequence vector, upper case,
    /// `CAlnVec::x_GetSeqVector`): residues of a protein, nucleotides of a nucleotide sequence,
    /// which `aln_vec::AlnMap::seq_string` reverse-complements or translates for a row.
    fn letters(&self, seq: SeqKey, from: u32, to: u32, out: &mut Vec<u8>);

    /// The genetic code of a translated row (`CAlnVec::SetGenCode`, showalign.cpp:2100-2101): the
    /// query's `m_MasterGeneticCode` for row 0 (`master`), the database's `m_SlaveGeneticCode` for
    /// every other row, whatever its sequence (a row of the query other than row 0 too: BLASTX's
    /// query rows of HSPs of length 0, TBLASTX's rows of a query record that is also a subject).
    /// `None`: a code NCBI's table does not have, BLASTX's database code 0
    /// (`CBlastOptionsLocal::SetProgram` clears it for a protein database), on which a
    /// translation throws (`aln_vec::UnknownGeneticCode`).
    fn genetic_code(
        &self,
        _master: bool,
    ) -> Option<crate::algorithm::tblastx::report::DisplayCode> {
        None
    }

    /// Appends the label of `seq`'s rows (`x_FillSeqid` without `eShowBlastStyleId`).
    fn label(&self, seq: SeqKey, out: &mut Vec<u8>);

    /// Appends the masks of `seq`: the query's `TMaskedQueryRegions`, a subject's of
    /// `GetSubjectMasks` (inclusive, 0-based, in the records), each with its frame (`Mask`).
    fn masks(&self, seq: SeqKey, out: &mut Vec<Mask>);
}

/// A mask (`CSeqLocInfo`): its first and last position (inclusive, 0-based, in the record) and its
/// frame (`GetFrame`). A translated query's masks have the frame of the context they were made in
/// (`BLAST_ContextToFrame`: 1 to 3, -1 to -3, counted in the searched letters,
/// `Blast_GetSeqLocInfoVector`, blast_aux.cpp:930-959); every other mask has none
/// (`CSeqLocInfo::eFrameNotSet`, 0: BLASTN's query masks, blast_aux.cpp:892; a protein query's
/// context frame 0; the subjects' masks, seqinfosrc_seqvec.cpp:142, seqinfosrc_seqdb.cpp:118).
pub(crate) type Mask = (u32, u32, i8);

/// What a program's report gives the sections of its queries: [`Sequences`] by the query's and
/// the subjects' indices in the report (a subject's masks are those of its HSP list with the
/// query, `GetSubjectMasks` of the query's results).
pub(crate) trait ReportSequences: Sync {
    /// Whether subject `s_idx` is query `q_idx`'s own Bioseq (`Hsps::set_query_bioseq`).
    fn subject_is_query(&self, q_idx: usize, s_idx: u32) -> bool;
    fn query_letters(&self, q_idx: usize, from: u32, to: u32, out: &mut Vec<u8>);
    fn subject_letters(&self, s_idx: u32, from: u32, to: u32, out: &mut Vec<u8>);
    fn query_label(&self, q_idx: usize, out: &mut Vec<u8>);
    fn subject_label(&self, s_idx: u32, out: &mut Vec<u8>);
    fn query_masks(&self, q_idx: usize, out: &mut Vec<Mask>);
    fn subject_masks(&self, q_idx: usize, s_idx: u32, out: &mut Vec<Mask>);

    /// The genetic code of the query's translated row 0 (`Sequences::genetic_code`).
    fn query_genetic_code(&self) -> Option<crate::algorithm::tblastx::report::DisplayCode> {
        None
    }

    /// The genetic code of the translated rows after row 0, the database's
    /// (`Sequences::genetic_code`); `None` for BLASTX's code 0.
    fn subject_genetic_code(&self) -> Option<crate::algorithm::tblastx::report::DisplayCode> {
        None
    }
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

    fn masks(&self, seq: SeqKey, out: &mut Vec<Mask>) {
        match seq {
            SeqKey::Query => self.sequences.query_masks(self.q_idx, out),
            SeqKey::Subject(s_idx) => self.sequences.subject_masks(self.q_idx, s_idx, out),
        }
    }

    fn genetic_code(&self, master: bool) -> Option<crate::algorithm::tblastx::report::DisplayCode> {
        if master {
            self.sequences.query_genetic_code()
        } else {
            self.sequences.subject_genetic_code()
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
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:914-930
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

/// Appends the letters `from..=to` (positions in the records) of a query whose letters are
/// `letters` (a record's in upper case, a PSSM query's Bioseq letters as they are:
/// `row_letters::QueryLetters::get`) and start at the record position `shift` (the searched
/// letters of TBLASTN and TBLASTX: the `-query_loc` interval's start, `Placements::offset`;
/// BLASTN's whole record: 0, `algorithm::blastn::pairwise::QueryRecords`): the letters that the
/// sections of the reports read (`ReportSequences::query_letters`), as far as the letters go.
pub(crate) fn query_record_letters(
    letters: crate::report::row_letters::QueryLetters<'_>,
    shift: usize,
    from: u32,
    to: u32,
    out: &mut Vec<u8>,
) {
    for position in from as usize - shift..=to as usize - shift {
        let Some(letter) = letters.get(position) else {
            break;
        };
        out.push(letter);
    }
}

/// The label of the rows of a FASTA record read as a query or a `-subject` subject: its Seq-ids
/// with `-parse_deflines` (`fill_seqid`), else its local id (`Query_N`, `Subject_N`).
pub(crate) fn record_label(
    record: &crate::blastinput::fasta_reader::FastaRecord,
    show_gis: bool,
    out: &mut Vec<u8>,
) {
    seq_id_label(record.believed(), record.seq_id(), show_gis, out);
}

/// `record_label` of a record whose `-parse_deflines` Seq-ids are `believed` (`None` without the
/// option) and whose Seq-id is `seq_id` (BLASTX's records, which have their own type).
pub(crate) fn seq_id_label(
    believed: Option<&crate::blastinput::fasta_reader::believed::BelievedIds>,
    seq_id: crate::report::seqid::SeqId<'_>,
    show_gis: bool,
    out: &mut Vec<u8>,
) {
    match believed {
        Some(ids) => fill_seqid(&ids.ids, show_gis, out),
        None => seq_id.sam_name(out),
    }
}

/// The label of the rows of a database sequence: its parsed Seq-ids
/// (`ParsedSubject::anchored_label`), else `gnl|BL_ORD_ID|<oid>`, whose label is the tag
/// (`GetLabel` of a `general` id, `fill_seqid`): its OID in the opened database, not its ordinal id
/// where a list or an alias range leaves sequences out (`seqdb::Oid`, as
/// `report::record::SubjectNames::Database`; review 3 of line RG T3a, D2).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2231-2236
pub(crate) fn database_label(
    parsed: Option<&crate::report::db_subject::ParsedSubject>,
    oid: crate::seqdb::Oid,
    out: &mut Vec<u8>,
) {
    match parsed {
        Some(parsed) => parsed.anchored_label(out),
        None => out.extend_from_slice(oid.get().to_string().as_bytes()),
    }
}

/// Whether the `-subject` record `subject` is query `query`'s own record: the query is one record
/// with a subject (the same definition line and letters, `BelievedIds::shared`), and the subject
/// has its id (any other subject with the query's id is refused before the reports).
pub(crate) fn subject_is_query_record(
    query: &crate::blastinput::fasta_reader::FastaRecord,
    subject: &crate::blastinput::fasta_reader::FastaRecord,
) -> bool {
    match (query.believed(), subject.believed()) {
        (Some(query), Some(subject)) => query.shared.is_some() && subject.best == query.best,
        _ => false,
    }
}

/// The refusal of a run whose query `record` (one record as a query and as a subject under
/// `-parse_deflines`) has a query-anchored section that NCBI BLAST+ aborts on (`ncbi_aborts`,
/// asked by `pairwise::section_aborts`), named as the refusal of its masks in outfmt 0 names it
/// (`believed::refuse_shared_masks`).
pub(crate) fn aborted_section(
    record: &crate::blastinput::fasta_reader::FastaRecord,
    program: &str,
) -> anyhow::Error {
    let mut fasta = Vec::new();
    if let Some(ids) = record.believed() {
        crate::seqid::label::fasta(&ids.best, &mut fasta);
    }
    anyhow::anyhow!(
        "the record {} is both a query and a subject under -parse_deflines and NCBI BLAST+ aborts \
         on its query-anchored alignments (outfmt 1 to 4), whose merge leaves out the subject's \
         rows; this is not supported by LOSAT's {program}",
        String::from_utf8_lossy(&fasta)
    )
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
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2055-2111
/// Each group's HSPs are added in order (`aln_mix::merge`; an HSP whose `Add` throws is left
/// out) and merged; `numDistinctFrames` counts the groups with HSPs, also one whose HSPs all
/// threw (none can: a sequence's first match never meets another strand). A group without a
/// merged Dense-seg is skipped (`GetDenseg` throws). A merge that throws (`Merge` is outside the
/// `try`: TBLASTX's query record that is also a subject, `aln_mix::MixError::NcbiThrows`) fails
/// the section with NCBI's exception before its first byte, as NCBI's leaves `DisplaySeqalign`
/// there (`NcbiException::into_io_error`).
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
    let groups = hsps.groups();
    let mut merged = Vec::new();
    for (group, indices) in groups.iter().enumerate() {
        if indices.is_empty() {
            continue;
        }
        let denseg = aln_mix::merge(hsps, indices).map_err(|error| match error {
            aln_mix::MixError::Internal(what) => internal(what),
            // `ncbi_aborts`, which the program's report asks before the query's first byte.
            aln_mix::MixError::NcbiAborts => internal("a section that NCBI BLAST+ aborts on"),
            // `section_frees`, asked as `ncbi_aborts` is.
            aln_mix::MixError::NcbiFreesSegment => {
                internal("a section whose merge NCBI BLAST+ frees a segment in")
            }
            // NCBI's exception, before the section's first byte (`Merge` comes before `out <<
            // "\n"`): the program's report ends with it (`cli::NativeError::of_writer`).
            aln_mix::MixError::NcbiThrows(exception) => exception.into_io_error(),
        })?;
        if let Some(denseg) = denseg {
            merged.push((group, denseg));
        }
    }
    let num_distinct_frames = groups.iter().filter(|indices| !indices.is_empty()).count();
    writer.write_all(b"\n")?;
    let mut buffers = aln_vec::MapBuffers::default();
    for (group, denseg) in &merged {
        let mut map = AlnMap::new_in(denseg, &mut buffers).map_err(internal)?;
        if num_distinct_frames > 1 {
            writer.write_all(b"For reading frame ")?;
            writer.write_all(K_FRAME_CONVERSION[*group].as_bytes())?;
            writer.write_all(b" of query sequence:\n\n")?;
        }
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1206-1208
        if view.anchored {
            map.set_anchor(0).map_err(internal)?;
        }
        // The query's row (row 0) has the frame of its group (`frame[row]`, `prepare_rows`).
        let query_frame = if hsps.query_width == 3 {
            K_GROUP_FRAMES[*group]
        } else {
            0
        };
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2092-2109
        // (cited at this function)
        // A row that `x_PrepareRowData` cannot translate (`aln_vec::UnknownGeneticCode`) throws
        // inside the `try` of the group, after its heading and before any line of it: NCBI
        // writes nothing more of the group.
        let rows = match prepare_rows(&map, line_length, hsps.query_bioseq, query_frame, sequences)
        {
            Err(error) if aln_vec::UnknownGeneticCode::is(&error) => continue,
            rows => rows?,
        };
        display_rows(writer, &map, view, line_length, rows, sequences)?;
    }
    Ok(())
}

/// The query's frame of each frame group (`k_FrameConversion` as numbers).
const K_GROUP_FRAMES: [i8; K_NUM_FRAME] = [1, 2, 3, -1, -2, -3];

/// Whether NCBI BLAST+ aborts on the section of `hsps` (`aln_mix::MixError::NcbiAborts`): its
/// merge leaves out the rows of the query's own Bioseq while segments hold them. A program asks
/// it before any format of the search's call is written and refuses the run there
/// (`pairwise::blastp_section_aborts`; in a block search, after the reports of the blocks
/// before): NCBI's output before the abort is what its stream had flushed, which depends on its
/// buffer (to a file or a pipe, the version line and the two `NcbiEndl` of `PrintProlog`,
/// blast_format.cpp:387-395, plus each full buffer; more on a terminal), so no prefix of it is
/// NCBI's and LOSAT writes none of the call's. Only a section with the query's own Bioseq can get
/// there (an ordinary subject's row takes its first match, `x_SecondRowFits`), so the others are
/// not merged here. A merge that throws first (`aln_mix::MixError::NcbiThrows`, TBLASTX's) is not
/// an abort: the application catches NCBI's exception and exits, its stream flushed, and
/// `write_section` writes what NCBI writes before it and fails with it.
pub(crate) fn ncbi_aborts(hsps: &Hsps) -> io::Result<bool> {
    if hsps.query_bioseq.is_none() || hsps.is_empty() {
        return Ok(false);
    }
    for indices in hsps.groups().iter().filter(|indices| !indices.is_empty()) {
        match aln_mix::merge(hsps, indices) {
            Ok(_) => {}
            Err(aln_mix::MixError::NcbiAborts) => return Ok(true),
            Err(aln_mix::MixError::NcbiThrows(_)) => return Ok(false),
            Err(aln_mix::MixError::NcbiFreesSegment) => {
                return Err(internal(
                    "a match of length 0 in a section with the query's own Bioseq",
                ))
            }
            Err(aln_mix::MixError::Internal(what)) => return Err(internal(what)),
        }
    }
    Ok(false)
}

/// Whether NCBI BLAST+'s merge of the section of `hsps` frees a segment that it still reads
/// (`aln_mix::MixError::NcbiFreesSegment`): with a BLASTX HSP of length 0
/// (`HspStrands::no_widths`, in the first group) or a BLASTP one (a segment of length 0), the merge's split of a segment of letters of
/// its query row at a length of 0 (`aln_mix::Mix::place`), or `FillUnalignedRegions` replacing
/// the last reference to a segment of length 0 that `m_Segments` lists
/// (`aln_mix::Mix::fill_unaligned_regions`). The groups are merged in order, all before the
/// section's first byte (`write_section`), so the first group's merge comes first; only a group
/// with such an HSP is merged here (only such HSPs make segments of length 0). The program asks
/// it before any byte of the call's reports, as `ncbi_aborts`, and refuses the report there
/// (PD-LOSAT-NCBI-DEFECTS rule 3): what NCBI writes before its crash is what its stream had
/// flushed.
pub(crate) fn section_frees(hsps: &Hsps) -> io::Result<bool> {
    for indices in hsps.groups().iter() {
        // Only an HSP with a segment of length 0 makes the merge free one: BLASTX's of length 0
        // (`no_widths`) and BLASTP's `-ungapped` one (`Hsps::push_dense_seg` of a Dense-seg of
        // one segment of length 0).
        if !indices.iter().any(|&index| {
            let hsp = hsps.hsp(index);
            hsp.strands.no_widths || hsp.segments.iter().any(|segment| segment.len == 0)
        }) {
            continue;
        }
        match aln_mix::merge(hsps, indices) {
            Ok(_) | Err(aln_mix::MixError::NcbiThrows(_)) => {}
            Err(aln_mix::MixError::NcbiFreesSegment) => return Ok(true),
            Err(aln_mix::MixError::NcbiAborts) => {
                return Err(internal(
                    "a section of length 0 matches that NCBI BLAST+ aborts on",
                ))
            }
            Err(aln_mix::MixError::Internal(what)) => return Err(internal(what)),
        }
    }
    Ok(false)
}

/// `SAlnRowInfo`: what `x_PrepareRowData` makes of each row.
struct RowData {
    sequence: Vec<Vec<u8>>,
    strings: Vec<aln_vec::RowString>,
    labels: Vec<Vec<u8>>,
    row_rng: Vec<(i64, i64)>,
    /// `masked_regions` and `frame`: each row's masks as alignment ranges with their frames, by
    /// the row's sequence (`masks`), and the frame of the row, whose masks it shows
    /// (`x_OutputSeq`); `None` for a subject's row of width 3, whose frame no mask of a subject
    /// has.
    masked: Vec<Option<(usize, i8)>>,
    masks: Vec<Vec<(i64, i64, i8)>>,
    max_id_len: usize,
    max_start_len: usize,
}

/// `x_PrepareRowData` for these formats: no features, no statistics, no sequence labels.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1242-1250
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1269-1283
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:411-422 (`s_GetFrame`)
/// A row of width 3 has a frame of 1 to 3 or -1 to -3, every other row 0. The query's row of a
/// translated query (BLASTX, TBLASTX) has the frame of its group (`query_frame`: its first and
/// last positions are those of the group's HSPs, all in that frame of the record), and shows the
/// query's masks of that frame. The subjects' masks have none (`Mask`), so `x_OutputSeq`
/// (`locFrame == frame`) shows none on a subject's row of width 3 (TBLASTN's, TBLASTX's).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1363-1370
/// (`maxStartLen` is the width of the largest 0-based position, which the 1-based numbers
/// printed can exceed by one digit.)
fn prepare_rows(
    map: &AlnMap<'_>,
    line_length: usize,
    query_bioseq: Option<u32>,
    query_frame: i8,
    sequences: &dyn Sequences,
) -> io::Result<RowData> {
    let num_rows = map.num_rows();
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
    // The index of each sequence's list in `data.masks`, made once for the section.
    let mut lists: std::collections::HashMap<SeqKey, usize> = std::collections::HashMap::new();
    let mut intervals = Vec::new();
    for row in 0..num_rows {
        let key = map.seq_key(row);
        data.row_rng.push(map.seq_aln_range(row));
        let mut string = map.whole_aln_seq_string(row, line_length as u32, sequences)?;
        data.sequence.push(std::mem::take(&mut string.buffer));
        data.strings.push(string);
        let mut label = Vec::new();
        sequences.label(key, &mut label);
        data.max_id_len = data.max_id_len.max(label.len());
        data.labels.push(label);
        let max_cood = map.seq_start(row).max(map.seq_stop(row));
        data.max_start_len = data.max_start_len.max(max_cood.to_string().len());
        // The rows of one sequence have one list (`x_FillLocList` of its masks), whose masks of
        // the row's frame a row shows.
        let index = *lists.entry(key).or_insert_with(|| {
            intervals.clear();
            sequences.masks(key, &mut intervals);
            data.masks
                .push(fill_loc_list(map, key, query_bioseq, &intervals));
            data.masks.len() - 1
        });
        let frame = if row == 0 {
            Some(query_frame)
        } else {
            (map.width(row) == 1).then_some(0)
        };
        data.masked.push(frame.map(|frame| (index, frame)));
    }
    Ok(data)
}

/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4291-4338
/// A mask goes to the columns of the first row of its sequence whose letters it meets (on a
/// minus-strand row its end is the first column), and `x_OutputSeq` puts those columns in lower
/// case in every row of that sequence of frame 0: a subject's mask spreads to the same columns of
/// its other HSPs' rows, also of the other strand (PD-LOSAT-NCBI-DEFECTS rule 2, DECISIONS
/// T3-6). The rows of a subject that is the query's own Bioseq (`query_bioseq`) have the query's
/// id, so the query's row, the first, takes its masks.
fn fill_loc_list(
    map: &AlnMap<'_>,
    key: SeqKey,
    query_bioseq: Option<u32>,
    masks: &[Mask],
) -> Vec<(i64, i64, i8)> {
    let own = |row_key: SeqKey| match query_bioseq {
        Some(subject) => {
            let one = |key: SeqKey| key == SeqKey::Query || key == SeqKey::Subject(subject);
            row_key == key || (one(row_key) && one(key))
        }
        None => row_key == key,
    };
    let mut list = Vec::new();
    for &(from, to, frame) in masks {
        for row in 0..map.num_rows() {
            if own(map.seq_key(row)) && map.seq_start(row) <= to && from <= map.seq_stop(row) {
                let (first, last) = if map.is_positive_strand(row) {
                    (from, to)
                } else {
                    (to, from)
                };
                list.push((
                    i64::from(map.aln_pos_from_seq_pos(row, first)),
                    i64::from(map.aln_pos_from_seq_pos(row, last)),
                    frame,
                ));
                break;
            }
        }
    }
    list
}

/// `x_DisplayRowData`: the lines of `line_length` columns.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1428-1468
/// (`j` stays below 2^31: an alignment of more than 2^30 columns would make NCBI's `int` wrap.
/// An alignment of no columns, the merge of BLASTX HSPs of length 0 alone, ends at `TSeqPos`
/// 4294967295, which `(int)` makes -1: no line.)
fn display_rows<W: Write + ?Sized>(
    writer: &mut W,
    map: &AlnMap<'_>,
    view: View,
    line_length: usize,
    mut rows: RowData,
    sequences: &dyn Sequences,
) -> io::Result<()> {
    let aln_stop = i64::from(map.aln_stop() as i32);
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
        rows.masked[row].map(|(index, frame)| (&rows.masks[index][..], frame)),
        false,
    );
    add_space(out, 2);
    if !empty {
        out.extend_from_slice(end.to_string().as_bytes());
    }
    out.push(b'\n');
}

/// `x_OutputSeq` of a row line: the columns `start..start + len` of `sequence`, those of the row's
/// masks in lower case (`m_SeqLocChar` is `eLowerCase`, blast_format.cpp:1581); under `html`
/// (`eHtml`) the columns of each mask on the line in a `<font>` span (`mask_spans`,
/// `write_styled`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2486-2545
/// (The row's masks are those of its own sequence with its frame, `prepare_rows`. The
/// mismatch colours of the first branch need `eColorDifferentBases`, which no BLAST program
/// sets: blast_format.cpp, survey T4 §1.1.)
fn output_seq(
    out: &mut Vec<u8>,
    sequence: &[u8],
    start: usize,
    len: usize,
    masks: Option<(&[(i64, i64, i8)], i8)>,
    html: bool,
) {
    let end = (start + len).min(sequence.len());
    let first = out.len();
    out.extend_from_slice(&sequence[start.min(end)..end]);
    let Some((masks, row_frame)) = masks else {
        return;
    };
    let mut spans = Vec::new();
    for &(from, to, _) in masks.iter().filter(|&&(_, _, frame)| frame == row_frame) {
        let lo = from.max(start as i64);
        let hi = to.min(start as i64 + len as i64 - 1);
        let mut column = lo;
        while column <= hi {
            if let Some(letter) = out.get_mut(first + (column - start as i64) as usize) {
                letter.make_ascii_lowercase();
            }
            column += 1;
        }
        if html {
            if let Some(span) = mask_span(from, to, start as i64, len as i64) {
                spans.push(span);
            }
        }
    }
    if !spans.is_empty() {
        let letters = out.split_off(first);
        write_styled(out, &letters, start as i64, &spans);
    }
}

/// `eachSeqloc` of a mask of the columns `from..=to` on the line of columns `start..start + len`
/// (`x_OutputSeq` under `eHtml`, cited at `output_seq`): its first column on the line, and its
/// last where the mask ends on the line, else 0 (the `GetTo()` of `CRange<int>(0, 0)`); `None` for
/// a mask that does not meet the line and for the range `(0, 0)` that NCBI drops: a mask that
/// starts at the alignment's first column and does not end on the first line (its letters are in
/// lower case, without a span), or that covers that column alone.
fn mask_span(from: i64, to: i64, start: i64, len: i64) -> Option<(i64, i64)> {
    let lo = from.max(start);
    let hi = to.min(start + len - 1);
    if lo > hi {
        return None;
    }
    // `i == min<int>(to, start+len)`: only the mask's own last column, `to`, when it lies on the
    // line (`hi` is below `start + len`).
    let last = if to <= hi { to } else { 0 };
    (!(lo == 0 && last == 0)).then_some((lo, last))
}

/// The letters `letters` of a line of first column `start` with the spans `spans` (`mask_span`)
/// each in NCBI's mask template: a span opens at its first column and closes after its last
/// column on the line, unless that column is the line's first (`to > 0`), in which case it runs,
/// as a span that ends past the line does, to the line's end (`s_ProcessStyledContent` closes it
/// at the last letter).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2541-2573
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:374-387
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:93,111,163
/// (The letters hold no `<@`, so the template's one replacement is the letters.)
fn write_styled(out: &mut Vec<u8>, letters: &[u8], start: i64, spans: &[(i64, i64)]) {
    const OPEN: &[u8] = b"<font color=\"#000000\">";
    const CLOSE: &[u8] = b"</font>";
    let (mut front, mut end) = (false, false);
    let mut styled: Vec<u8> = Vec::new();
    for (i, &letter) in letters.iter().enumerate() {
        let i = i as i64;
        for &(from, to) in spans {
            if from - start == i {
                front = true;
            }
            if to - start == i && to - start > 0 {
                end = true;
            }
        }
        let (start_styled, stop_styled) = (front, end && front);
        let mut is_styled = false;
        if start_styled || (!styled.is_empty() && !stop_styled) {
            styled.push(letter);
            is_styled = true;
        }
        if !styled.is_empty() && (stop_styled || i == letters.len() as i64 - 1) {
            out.extend_from_slice(OPEN);
            out.append(&mut styled);
            out.extend_from_slice(CLOSE);
        }
        if !is_styled {
            out.push(letter);
        }
        if end && front {
            end = false;
            front = false;
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
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3099-3120
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3076-3096
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
        if let Some(byte) = column(&mut line, insert.aln_start - aln_start + 1)? {
            *byte = b'\\';
        }
    }
    let mut inserts = Vec::new();
    do_fills(
        map,
        row,
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

/// The byte of `line` at `index` where NCBI writes the mark of an insertion (`\` or `|`); `None`
/// at `index == line.len()`, where NCBI's write changes no byte of the line.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3089 (`x_FillInserts`)
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3050 and 3060
/// (`x_DoFills`, both branches)
/// Both strings hold `aln_range.GetLength()` bytes, and `x_GetInserts` takes an insertion whose
/// alignment start is below the line's end, so the index is at most the length. It is the length
/// for an insertion after the last column of the alignment's last line: a row that ends with an
/// insertion, which the merge puts after the alignment's last column (the split row of the query's
/// own Bioseq, `aln_mix::Mix::place`, whose HSP's longest segment follows a subject insertion).
/// NCBI then writes `\` and `|` onto the strings' terminator (`operator[]` at `size()`): no byte
/// of the line changes, the lines print their spaces, and `s_AdjustInsert` puts the letters past
/// the line, right-aligned at that column. That visible result (no mark) is deterministic in NCBI
/// BLAST+ 2.17.0 (review 2 of line RG T3a: rows `db.of1.bp.lo24.f1`, `lp17.f2`, `lp58.f2`, and
/// `lp20.ll20`, where `-line_length 20` puts the insertion on no line) and LOSAT reproduces it
/// (PD-LOSAT-NCBI-DEFECTS rule 2). An index past the length, which `x_GetInserts` cannot give, is
/// LOSAT's internal error.
fn column(line: &mut [u8], index: i64) -> io::Result<Option<&mut u8>> {
    let len = line.len();
    match usize::try_from(index) {
        Ok(index) if index == len => Ok(None),
        Ok(index) if index < len => Ok(Some(&mut line[index])),
        _ => Err(internal("an insertion outside its line")),
    }
}

/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3031-3073
/// (The recursion as a loop.) The letters of an insertion are not masked; they are the row's
/// (`GetSeqString`: reverse-complemented on a minus-strand row, translated on a row of width 3,
/// whose insertion's length is in nucleotides).
fn do_fills(
    map: &AlnMap<'_>,
    row: usize,
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
            if let Some(byte) = column(&mut bar, cur_insert_aln_start - aln_start + 1)? {
                *byte = b'|';
            }
            if is_first_insert || cur_insert_aln_start - prvs_insert_aln_end >= 1 {
                new_insert.clear();
                let seq_end = insert.seq_start + insert.insert_len - 1;
                map.seq_string(row, insert.seq_start, seq_end, sequences, &mut new_insert)?;
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
#[cfg(test)]
mod tests_t3b;
#[cfg(test)]
mod tests_t3c;
#[cfg(test)]
mod tests_t3c_self;
#[cfg(test)]
mod tests_t4a0;
#[cfg(test)]
mod tests_x2g;
