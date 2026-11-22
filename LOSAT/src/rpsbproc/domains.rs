//! The domain hits of one protein query: the alignments of an archive turned into rpsbproc's
//! `TDomSeqAlignment`s, their classification into regions, representative, standard and full
//! sets (`CCdAlignProcessor::Calculate`), and the records that the DOMAINS output prints
//! (`TDomSeqAlignIndex::CreateRecordSets`).

use std::cmp::Ordering;
use std::collections::HashMap;

use super::archive::{Align, DenseSeg, StdLoc, StdPair};
use super::data::{dom_acc_type, sort_ties_in_input_order, CddData, TOTALSIGS};
use super::frames::{minus_pr_to_plus_na, plus_pr_to_plus_na, plus_rf_idx_to_minus, RF_SIZE};
use super::segset::SegSet;
use super::sites::ResiduePos;

/// `TDomSeqAlignment::GAP_THRESHOLD`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:1203-1205
const GAP_THRESHOLD: i64 = 35;

/// `TDomSeqAlignment`, with the members that the domain output reads.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:1084-1123,1206-1221
// (`m_iFrom` and `m_iTo` are `i64`: the XML input's `GetQuery_from` less `COORDSBASE`, for a value
// of `-2147483648`, stores one less than `INT_MIN` there, the value of NCBI's overflowing `int`
// arithmetic, PD-LOSAT-NCBI-DEFECTS exception 12; an archive's are `int`s, converted from
// `TSeqPos`.)
#[derive(Clone, Debug)]
pub(crate) struct DomAlign {
    pub(crate) pssm_id: u32,
    pub(crate) from: i64,
    pub(crate) to: i64,
    pub(crate) m_starts: Vec<i32>,
    pub(crate) s_starts: SubjectStarts,
    pub(crate) lens: Vec<u32>,
    pub(crate) e_value: f64,
    pub(crate) bit_score: f64,
    pub(crate) region: i32,
    pub(crate) n_missing: f64,
    pub(crate) c_missing: f64,
    pub(crate) spec_qualified: bool,
    pub(crate) rep_class: i32,
    pub(crate) rep: bool,
    /// `m_ClipSet`: the query ranges (as positions on the model, `ori_from` side) that a
    /// structural-domain hit keeps for its motifs (set by `calculate`, empty otherwise).
    pub(crate) clip_set: SegSet,
    /// `m_eAlignType` being `ePr2Na`: a hit of a translated query (a `Std-seg` whose query locations
    /// carry a strand); its `m_starts` are positions on a "pseudo-protein".
    pub(crate) pr2na: bool,
    /// `m_ReadingFrame` (`TFRAMEINDEX`, an `unsigned int`): the plus-strand frame index of
    /// `from` in its two low bits, and for a minus-strand hit the aligned length in codons
    /// shifted left by 2.
    pub(crate) reading_frame: u32,
}

impl DomAlign {
    /// A hit without segments (`TSeqAlignment::Reset` and the start of `ParseAlignSegs`).
    fn empty(align: &Align, pssm_id: u32) -> Self {
        Self::new(
            pssm_id,
            i64::from(i32::MAX),
            0,
            align.e_value.unwrap_or(0.0),
            align.bit_score.unwrap_or(0.0),
        )
    }

    /// A hit of the range `from`..`to` without segments, of a protein query (the alignment
    /// type and the reading frame of `TDomSeqAlignment`'s constructor).
    pub(crate) fn new(pssm_id: u32, from: i64, to: i64, e_value: f64, bit_score: f64) -> Self {
        DomAlign {
            pssm_id,
            from,
            to,
            m_starts: Vec::new(),
            s_starts: SubjectStarts::default(),
            lens: Vec::new(),
            e_value,
            bit_score,
            region: 0,
            n_missing: 0.0,
            c_missing: 0.0,
            spec_qualified: false,
            rep_class: 0,
            rep: false,
            clip_set: SegSet::new(),
            pr2na: false,
            reading_frame: 0,
        }
    }

    /// The aligned length in codons of a minus-strand hit (`m_ReadingFrame` shifted right by 2), 0
    /// for the others.
    pub(crate) fn aligned_len(&self) -> u32 {
        self.reading_frame >> 2
    }

    /// The plus-strand frame index (`m_ReadingFrame` masked with `READINGFRAME::RF_SIZE`).
    pub(crate) fn plus_rf_idx(&self) -> u32 {
        self.reading_frame & RF_SIZE
    }
}

/// `m_vecSStarts` (`std::vector<SeqPos_t>`): the starts of the segments on the model. NCBI keeps
/// `int`s; LOSAT keeps `i32`s, 4 bytes a segment, until a start past `int`'s range comes, and
/// `i64`s from then on: `ParseAlignStrings` of the XML input moves its `int` `hf` past the range
/// with an increment and a decrement of `hf` (and `COORDSBASE` taken off an `Hsp_hit-from` of
/// `-2147483648`), which LOSAT computes as their values (PD-LOSAT-NCBI-DEFECTS exception 12). An
/// archive's starts are `int`s (converted from `TSeqPos`), so they stay narrow, and so does every
/// start of a real document.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:1121
// (No equality of two of them: the storage is not part of the value. The tests compare the values
// with a list, `PartialEq<Vec<i64>>`.)
#[derive(Clone, Debug)]
pub(crate) enum SubjectStarts {
    Narrow(Vec<i32>),
    Wide(Vec<i64>),
}

impl Default for SubjectStarts {
    fn default() -> Self {
        SubjectStarts::Narrow(Vec::new())
    }
}

impl SubjectStarts {
    /// `push_back`.
    pub(crate) fn push(&mut self, start: i64) {
        match self {
            SubjectStarts::Narrow(starts) => match i32::try_from(start) {
                Ok(start) => starts.push(start),
                Err(_) => {
                    let mut wide: Vec<i64> = starts.iter().map(|&start| i64::from(start)).collect();
                    wide.push(start);
                    *self = SubjectStarts::Wide(wide);
                }
            },
            SubjectStarts::Wide(starts) => starts.push(start),
        }
    }

    pub(crate) fn len(&self) -> usize {
        match self {
            SubjectStarts::Narrow(starts) => starts.len(),
            SubjectStarts::Wide(starts) => starts.len(),
        }
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// `m_vecSStarts[i]` (`i` below `len()`).
    pub(crate) fn get(&self, i: usize) -> i64 {
        match self {
            SubjectStarts::Narrow(starts) => i64::from(starts[i]),
            SubjectStarts::Wide(starts) => starts[i],
        }
    }
}

/// The starts' values against a list, whatever the storage (for the tests).
#[cfg(test)]
impl PartialEq<Vec<i64>> for SubjectStarts {
    fn eq(&self, other: &Vec<i64>) -> bool {
        self.len() == other.len()
            && other
                .iter()
                .enumerate()
                .all(|(i, &start)| self.get(i) == start)
    }
}

/// `TDomSeqAlignIndex`: the indices of the alignments in each set.
#[derive(Clone, Debug, Default)]
pub(crate) struct AlignIndex {
    pub(crate) sorted: Vec<usize>,
    pub(crate) concise: Vec<usize>,
    pub(crate) standard: Vec<usize>,
    pub(crate) qualified_feat: Vec<usize>,
    /// `m_vecSDIndice`: the structural-domain hits that keep a part of the query.
    pub(crate) structural: Vec<usize>,
}

/// `TDataModes`: `rep`, `std`, `full`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum DataMode {
    Rep,
    Std,
    Full,
}

/// `TDomSeqAlignIndex::__TCdAlignRecord` of the DOMAINS output.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Record {
    /// The alignment (`pAlign`).
    pub(crate) align: usize,
    /// The domain (`pCdInfo`, index in `CddData::domains`).
    pub(crate) domain: usize,
    /// The superfamily (`pClst`, index in `CddData::clusters`).
    pub(crate) cluster: Option<usize>,
}

/// `CleanAlignment`: the segments where both rows are aligned.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/basedata_impl.hpp:2-19
// (`starts` shorter than twice `lens` is read past its end by NCBI: refused by the caller.)
fn clean_alignment(seg: &DenseSeg, align: &mut DomAlign) {
    for (i, &len) in seg.lens.iter().enumerate() {
        let (m, s) = (seg.starts[2 * i], seg.starts[2 * i + 1]);
        if m >= 0 && s >= 0 {
            align.m_starts.push(m);
            align.s_starts.push(i64::from(s));
            align.lens.push(len);
        }
    }
}

/// `ParseAlignSegs` of a `Dense-seg` and `ParseAlignScores`: the alignment, without its
/// PSSM-ID (`ProcessCDAlign` sets it).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:9-43
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:214-262
// (`ParseAlignScores`, cited at `archive::ArchiveStream::align_score`, which keeps the scores
// it reads.)
// (One less than the sum of an entry of `m_vecMStarts` (an `int`) and of `m_vecLens` (an
// `unsigned int`) is computed in `unsigned int`: the sum wraps
// as an `unsigned int` and is stored back in an `int`.)
pub(crate) fn parse_align(seg: &DenseSeg, align: &Align, pssm_id: u32) -> DomAlign {
    let mut dom = DomAlign::empty(align, pssm_id);
    clean_alignment(seg, &mut dom);
    if let (Some(&first), Some(&last), Some(&last_len)) =
        (dom.m_starts.first(), dom.m_starts.last(), dom.lens.last())
    {
        dom.from = i64::from(first);
        dom.to = i64::from((last as u32).wrapping_add(last_len).wrapping_sub(1) as i32);
    }
    dom
}

/// `ParseAlignSegs` of the `Std-seg`s of a translated search and `ParseAlignScores`: the
/// alignment, without its PSSM-ID (`ProcessCDAlign` sets it), or the refusal of what NCBI
/// reads outside its vectors.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:44-133
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:166-209
// The types decide the arithmetic: `GetFrom`/`GetTo`/`GetPoint` are `TSeqPos` (unsigned) stored in
// `SeqPos_t` (`int`, wraps in NCBI's build); `t - f + 1` is `int` arithmetic, which overflows for
// an interval of 2^31 positions or more (undefined): computed as its value (PD-LOSAT-NCBI-DEFECTS
// exception 12), which `push_back` converts to the `unsigned int` of `m_vecLens` (modulo 2^32);
// `m_iFrom` modulo `RF_SIZE`, `*iterMS` divided by `RF_SIZE` and the minus-strand sum are `unsigned
// int` arithmetic (wraps); `m_iTo` plus 1 is `int` (`m_iTo` is at least 0), which overflows for a
// hit that ends at 2^31 - 1: computed as its value (exception 12), which the division by the
// `unsigned int` `RF_SIZE` converts modulo 2^32. A segment whose subject location is neither `int`
// nor `pnt` leaves its query start without a subject start and length: the vectors differ in
// length, which the plus-strand code and the classification read within bounds (LOSAT keeps the
// three vectors as NCBI does), but the minus-strand loop above reads `m_vecLens` past its end
// (refused). `ids` with fewer than 2 ids in the first segment is read past its end (refused, by the
// caller). The test of `seqlocs[0]` and `seqlocs[1]` for `null` (`IsNull`, joined by "or") reads
// `seqlocs[1]` only when the first location is not `null`: a segment without locations, or whose only location is not `null`, is
// read past its end (refused), and a segment whose first location is `null` is skipped whatever
// follows (PS16c review oracle `nullone`). The strand that counts is the last one set; `unknown`,
// `plus` and `both` make a plus-strand hit, `minus` and `both-rev` a minus-strand one, `other` (or
// none) an alignment of two proteins.
pub(crate) fn parse_std_align(
    segs: &[StdPair],
    align: &Align,
    pssm_id: u32,
) -> Result<DomAlign, String> {
    const STRAND_MINUS: i32 = 2;
    const STRAND_BOTH_REV: i32 = 4;
    const STRAND_OTHER: i32 = 255;
    let mut dom = DomAlign::empty(align, pssm_id);
    let mut qstrand = STRAND_OTHER;
    for &StdPair { query, hit } in segs {
        if query == StdLoc::Missing {
            return Err(
                "a Std-seg without locations (NCBI reads past them) is not supported by LOSAT's rpsbproc".to_string(),
            );
        }
        if query == StdLoc::Null {
            continue;
        }
        if hit == StdLoc::Missing {
            return Err(
                "a Std-seg whose only location is not null (NCBI reads past it) is not supported by LOSAT's rpsbproc".to_string(),
            );
        }
        if hit == StdLoc::Null {
            continue;
        }
        if query == StdLoc::Empty || hit == StdLoc::Empty {
            continue;
        }
        let (start, end, strand) = match query {
            StdLoc::Int { from, to, strand } => (from as i32, to as i32, strand),
            StdLoc::Pnt { point, strand } => (point as i32, point as i32, strand),
            _ => continue,
        };
        if i64::from(start) < dom.from {
            dom.from = i64::from(start);
        }
        if i64::from(end) > dom.to {
            dom.to = i64::from(end);
        }
        dom.m_starts.push(start);
        if let Some(strand) = strand {
            qstrand = strand;
        }
        match hit {
            StdLoc::Int { from, to, .. } => {
                let (f, t) = (i64::from(from as i32), i64::from(to as i32));
                dom.s_starts.push(f);
                dom.lens.push((t - f + 1) as u32);
            }
            StdLoc::Pnt { point, .. } => {
                dom.s_starts.push(i64::from(point as i32));
                dom.lens.push(1);
            }
            _ => {}
        }
    }
    if qstrand != STRAND_OTHER {
        dom.pr2na = true;
        dom.reading_frame = (dom.from as u32) % RF_SIZE;
        if qstrand == STRAND_MINUS || qstrand == STRAND_BOTH_REV {
            let aligned_len = (dom.to + 1) as u32 / RF_SIZE;
            dom.reading_frame |= aligned_len << 2;
            if dom.m_starts.len() > dom.lens.len() {
                return Err(
                    "a minus-strand Std-seg with a subject location other than int or pnt (NCBI reads past its segment lengths) is not supported by LOSAT's rpsbproc".to_string(),
                );
            }
            for (start, &len) in dom.m_starts.iter_mut().zip(&dom.lens) {
                *start =
                    aligned_len.wrapping_sub((*start as u32 / RF_SIZE).wrapping_add(len)) as i32;
            }
        } else {
            for start in &mut dom.m_starts {
                *start = (*start as u32 / RF_SIZE) as i32;
            }
        }
    }
    Ok(dom)
}

/// `TSeqAlignment::CalcMasterGaps` with `GAP_THRESHOLD`: the query ranges between segments that are
/// at least `GAP_THRESHOLD` long.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1324-1339
// (`m_vecLens[i]` cast to `SeqPos_t` converts an `unsigned int` to `int`, which wraps in NCBI's build;
// the `int` sums after it are computed in `i64`, `segset::Seg`.)
fn master_gaps(align: &DomAlign) -> SegSet {
    let mut gaps = SegSet::new();
    for i in 0..align.lens.len().saturating_sub(1) {
        let gap_start = i64::from(align.m_starts[i]) + i64::from(align.lens[i] as i32);
        let gap_len = i64::from(align.m_starts[i + 1]) - gap_start;
        if gap_len >= GAP_THRESHOLD {
            gaps.add(gap_start, gap_start + gap_len - 1);
        }
    }
    gaps
}

/// `TDomSrcCount`: how many hits of each source database a region took (one each).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:115,154-166
fn count_source(counts: &mut [i32; 7], accession: &[u8]) -> bool {
    const E_PRK: usize = 6;
    let source = dom_acc_type(accession);
    if source == TOTALSIGS {
        return true;
    }
    let source = source.min(E_PRK);
    if counts[source] >= 1 {
        return false;
    }
    counts[source] += 1;
    true
}

/// A hit as `Calculate` sorts it (`TDomAlignFacility`).
#[derive(Clone, Copy, Debug)]
struct Facility {
    align: usize,
    domain: usize,
    /// The gaps (`ulGapsIdx`, an index in the gap sets).
    gaps: Option<usize>,
}

/// `CCdAlignProcessor::Calculate` of the alignments `selected` (all of a protein query's, or
/// those of one reading frame of a nucleotide query's): their scores against their models, the
/// regions, and the sets of `index`. Structural domains (`sd`) are set apart and take no region
/// (their clipping, `m_vecSDIndice` and `m_ClipSet`, serves the motif output).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:48-130
// The models missing from the data are looked up again after `x_LoadMissingDomains`, which the
// offline processor does not override (`cdalignproc_base.hpp:123`, an empty body): they stay
// missing and their hits are dropped (cdalignproc_base.cpp:133-195). `Multidom` hits (`m_iRepClass`
// 1) are not made since 2018-04-04 (the commented lines above), so the output's `Multidom` type is
// never printed. (`m_uiLength` less the sum of `sstart` and `len` is `unsigned int` arithmetic
// stored in an `int`, of the start's value modulo 2^32 where the XML input took it past `int`
// (`SubjectStarts`); `iNMissing` is that value (PD-LOSAT-NCBI-DEFECTS exception 12); the ratios
// divide doubles, so a model of length 0 gives an infinity or a NaN as NCBI's does. The two sorts
// by `TDomAlignOrderByEValue` are Rust's standard sort with that order made total (`e_value_order`:
// a NaN E-value or bit score, which NCBI's comparison orders against no number, sorts after every
// number) and the hit's index as the last key (`sort_ties_in_input_order`; `selected` is in
// increasing order): hits equal in E-value and bit score keep the order of `selected`, where NCBI's
// `std::sort` may leave them in another, and so may become the representative or `Specific` hit in
// NCBI's place (PD-LOSAT-NCBI-DEFECTS 1.67, the order of tied hits).)
pub(crate) fn calculate(aligns: &mut [DomAlign], selected: &[usize], data: &CddData) -> AlignIndex {
    let mut index = AlignIndex::default();
    let mut gap_sets: Vec<SegSet> = Vec::new();
    let mut non_multi: Vec<Facility> = Vec::new();
    let mut structural: Vec<Facility> = Vec::new();
    for &i in selected {
        let align = &mut aligns[i];
        if align.s_starts.is_empty() {
            continue;
        }
        let Some(domain_index) = data.find_domain(align.pssm_id) else {
            continue;
        };
        let domain = &data.domains[domain_index];
        let length = domain.length;
        let n_missing = align.s_starts.get(0);
        let c_missing = length.wrapping_sub(
            (align.s_starts.get(align.s_starts.len() - 1) as u32)
                .wrapping_add(*align.lens.last().expect("not empty")),
        ) as i32;
        align.n_missing = n_missing as f64 / f64::from(length);
        align.c_missing = f64::from(c_missing) / f64::from(length);
        align.spec_qualified =
            domain.min_bit_score > 0.0 && align.bit_score >= domain.min_bit_score;
        let gaps = master_gaps(align);
        let gaps = if gaps.is_empty() {
            None
        } else {
            gap_sets.push(gaps);
            Some(gap_sets.len() - 1)
        };
        let facility = Facility {
            align: i,
            domain: domain_index,
            gaps,
        };
        if domain.is_struct_dom {
            align.rep_class = 2;
            structural.push(facility);
        } else {
            align.rep_class = 0;
            non_multi.push(facility);
        }
    }
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.hpp:62-70
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:197-270
    sort_ties_in_input_order(
        &mut non_multi,
        |facility| facility.align,
        |a, b| e_value_order(&aligns[a.align], &aligns[b.align]),
    );
    let segs_of = |facility: &Facility, aligns: &[DomAlign], gap_sets: &[SegSet]| {
        let mut segs = SegSet::new();
        segs.add(aligns[facility.align].from, aligns[facility.align].to);
        if let Some(gaps) = facility.gaps {
            segs.clip(&gap_sets[gaps]);
        }
        segs
    };
    let mut concise: Vec<Facility> = Vec::new();
    let mut source_counts: HashMap<i32, [i32; 7]> = HashMap::new();
    for facility in &non_multi {
        let this_segs = segs_of(facility, aligns, &gap_sets);
        let this_length = this_segs.total_residues() as f64;
        let mut overlaps = SegSet::new();
        let mut joined = false;
        for rep in &concise {
            let rep_segs = segs_of(rep, aligns, &gap_sets);
            let rep_length = rep_segs.total_residues() as f64;
            let mut overlap = this_segs.clone();
            overlap.cross(&rep_segs);
            let overlap_length = overlap.total_residues() as f64;
            if overlap_length > 0.0
                && (overlap_length / this_length > 0.5
                    || (overlap_length / rep_length > 0.5
                        && data.domains[facility.domain].cluster_id
                            == data.domains[rep.domain].cluster_id))
            {
                let region = aligns[rep.align].region;
                let align = &mut aligns[facility.align];
                align.region = region;
                align.rep = false;
                joined = true;
                break;
            }
            overlaps.merge(&overlap);
        }
        if !joined {
            let rep = overlaps.total_residues() as f64 / this_length < 0.5;
            aligns[facility.align].rep = rep;
            if rep {
                aligns[facility.align].region = concise.len() as i32;
                concise.push(*facility);
            }
        }
        let domain = &data.domains[facility.domain];
        let counts = source_counts
            .entry(aligns[facility.align].region)
            .or_insert([0; 7]);
        if count_source(counts, data.text(domain.accession)) {
            index.standard.push(facility.align);
            aligns[facility.align].rep = true;
            if domain.curated {
                index.qualified_feat.push(facility.align);
            }
        }
        index.sorted.push(facility.align);
    }
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:275-303
    // (`vecPrivileged` is always empty.)
    index.concise = concise.iter().map(|facility| facility.align).collect();
    if !index.sorted.is_empty() && index.standard.is_empty() {
        index.standard = index.sorted.clone();
    }
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:305-342
    // (`m_uiLength` less 1 is `unsigned int`, converted to `SeqPos_t` by `AddSeg`: a model of length
    // 0 gives -1, an empty range.)
    if !structural.is_empty() {
        sort_ties_in_input_order(
            &mut structural,
            |facility| facility.align,
            |a, b| e_value_order(&aligns[a.align], &aligns[b.align]),
        );
        let mut covered = SegSet::new();
        for facility in &structural {
            let mut sd_segs = SegSet::new();
            sd_segs.add(
                0,
                i64::from(data.domains[facility.domain].length.wrapping_sub(1) as i32),
            );
            map_seg_set(&mut sd_segs, &aligns[facility.align]);
            sd_segs.clip(&covered);
            if !sd_segs.is_empty() {
                let mut clip = SegSet::new();
                for seg in sd_segs.segs() {
                    clip.add(seg.ori_from, SegSet::ori_to(seg));
                }
                aligns[facility.align].clip_set = clip;
                covered.merge(&sd_segs);
                index.structural.push(facility.align);
            }
        }
    }
    index
}

/// `TDomAlignOrderByEValue` made a total order: the E-value ascending, then the bit score
/// descending, the numbers compared as C++'s `<` and `>` compare doubles (-0.0 and 0.0 are equal);
/// a NaN, which NCBI's comparison orders against no number, after every number, and two NaNs
/// equal. Hits equal in both keep their order (`sort_ties_in_input_order`); NCBI's `std::sort` may
/// leave them in another (PD-LOSAT-NCBI-DEFECTS 1.67, the order of tied hits).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.hpp:62-70
fn e_value_order(a: &DomAlign, b: &DomAlign) -> Ordering {
    key_order(a.e_value, b.e_value, false).then_with(|| key_order(a.bit_score, b.bit_score, true))
}

/// One key of `e_value_order`: two numbers in C++'s order (`descending`: reversed), a NaN after
/// every number.
fn key_order(x: f64, y: f64, descending: bool) -> Ordering {
    match (x.is_nan(), y.is_nan()) {
        (false, false) => {
            let order = x.partial_cmp(&y).expect("two numbers are ordered");
            if descending {
                order.reverse()
            } else {
                order
            }
        }
        (x_nan, y_nan) => x_nan.cmp(&y_nan),
    }
}

/// `CCdAlignProcessor::SortReadingFrames`: the indices of the alignments of each reading frame
/// (`GetRFIdx` of the query's length), or NCBI's exception for a query of length 0.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:20-45
pub(crate) fn sort_reading_frames(
    aligns: &[DomAlign],
    query_len: u32,
) -> Result<[Vec<usize>; 6], &'static str> {
    if query_len == 0 {
        return Err("Cannot sort reading frames - query sequence length unknown");
    }
    let mut frames: [Vec<usize>; 6] = Default::default();
    for (i, align) in aligns.iter().enumerate() {
        frames[rf_index(align, query_len) as usize].push(i);
    }
    Ok(frames)
}

/// `TSeqAlignment::GetRFIdx` of `seqLen`: the reading frame index (0 to 5) of an alignment.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1342-1351
pub(crate) fn rf_index(align: &DomAlign, query_len: u32) -> u32 {
    if !align.pr2na {
        return 0;
    }
    let rfidx = align.plus_rf_idx();
    if align.aligned_len() > 0 {
        plus_rf_idx_to_minus(rfidx, query_len)
    } else {
        rfidx
    }
}

/// `TSeqAlignment::MapSegSet` with `segset`, for a protein alignment.
pub(crate) fn map_seg_set(set: &mut SegSet, align: &DomAlign) {
    set.map_through(&align.m_starts, &align.s_starts, &align.lens);
}

/// `TSeqAlignment::GetTranslatedPosMap`: the residues of a mapped set, each with its position
/// on the query (on the translation of its frame for a hit of a translated query) and its
/// original position; `Err` with the text of NCBI's exception when a minus-strand hit reaches
/// past the query's length.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1252-1294
// (`qLen` less `rfidx` and `c` plus `offset` are `unsigned int` arithmetic, which wraps; the sum
// of `c` and `offset` is stored in an `int`. The increment of `c` and the sum of `ori_from` and
// `inc` are `int` arithmetic, which overflows past `INT_MAX` (undefined; a mapped range at the limit, made by
// hand): computed as their values, which the plus strand stores as they are (PD-LOSAT-NCBI-DEFECTS
// exception 12); a position outside the translation is then read by the out-of-range rules of its
// reader.)
pub(crate) fn translated_pos_map(
    set: &SegSet,
    align: &DomAlign,
    query_len: u32,
) -> Result<Vec<ResiduePos>, String> {
    // `Some(offset)` on the minus strand, whose `c + offset` is an `unsigned int`.
    let mut minus_offset = None;
    if align.pr2na {
        let aligned_len = align.aligned_len();
        let aa_len = query_len.wrapping_sub(align.plus_rf_idx()) / RF_SIZE;
        if aligned_len > 0 {
            if aligned_len > aa_len {
                return Err(format!(
                    "Invalid protein length {aa_len}: shorter than aligned range {aligned_len}"
                ));
            }
            minus_offset = Some(aa_len - aligned_len);
        }
    }
    let mut positions = Vec::with_capacity(set.total_residues().max(0) as usize);
    for seg in set.segs() {
        let mut c = seg.from;
        let mut inc = 0i64;
        while c <= seg.to {
            let curr = match minus_offset {
                Some(offset) => i64::from((c as u32).wrapping_add(offset) as i32),
                None => c,
            };
            positions.push(ResiduePos {
                curr,
                ori: seg.ori_from + inc,
            });
            c += 1;
            inc += 1;
        }
    }
    Ok(positions)
}

/// The two ends of a mapped set after `TSeqAlignment::Pr2NaConvert` (`GetLeft`, `GetRight`):
/// positions on the query's plus strand for a hit of a translated query; the set's own ends
/// for the others (`None`: an empty set). A position that NCBI's `int` arithmetic took past
/// `int`'s range holds its value, and `pr - 1` of `INT_MIN` is computed as its value
/// (`frames::minus_pr_to_plus_na`; PD-LOSAT-NCBI-DEFECTS exception 12): the conversions after
/// them are `unsigned int`'s, modulo 2^32.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1137-1177
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/segset_impl.hpp:281-290
// Only the ends are printed (the motif lines), so only the first and the last ranges are converted;
// `ori_from` is not read after the conversion. The minus strand reverses the list: the left end
// comes from the last range's `to`, the right end from the first range's `from`. The minus strand's
// update of `ori_from` (it grows by the range length divided by `m_iFactor`, less one) is not
// computed: its value is never read, so no printed number depends on it (PS16c review S3-4).
pub(crate) fn pr2na_ends(set: &SegSet, align: &DomAlign) -> Option<(i64, i64)> {
    let (first, last) = (set.segs().first()?, set.segs().last()?);
    if !align.pr2na {
        return Some((first.from, last.to));
    }
    let rfidx = align.plus_rf_idx();
    let aligned_len = align.aligned_len();
    if aligned_len > 0 {
        let na_len = aligned_len.wrapping_mul(RF_SIZE).wrapping_add(rfidx);
        let left = (minus_pr_to_plus_na(last.to, RF_SIZE, na_len) as u32)
            .wrapping_sub(RF_SIZE)
            .wrapping_add(1) as i32;
        let right = minus_pr_to_plus_na(first.from, RF_SIZE, na_len);
        Some((i64::from(left), i64::from(right)))
    } else {
        let left = plus_pr_to_plus_na(first.from, rfidx);
        let right = (plus_pr_to_plus_na(last.to, rfidx) as u32)
            .wrapping_add(RF_SIZE)
            .wrapping_sub(1) as i32;
        Some((i64::from(left), i64::from(right)))
    }
}

/// `TDomSeqAlignIndex::CreateRecordSets`, the domain records (`rDomAligns`) of `mode`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1918-2010
// The hierarchy roots and the site records (`rFeatAligns`, `ExtractFeatAligns`) serve the site
// output (`feature_records`).
pub(crate) fn domain_records(
    index: &AlignIndex,
    aligns: &[DomAlign],
    data: &CddData,
    mode: DataMode,
) -> Vec<Record> {
    let indices = match mode {
        DataMode::Rep => &index.concise,
        DataMode::Std => &index.standard,
        DataMode::Full => &index.sorted,
    };
    indices
        .iter()
        .map(|&align| {
            let domain = data
                .find_domain(aligns[align].pssm_id)
                .expect("a classified hit has its model");
            let cluster_id = data.domains[domain].cluster_id;
            let cluster = if cluster_id > 0 {
                data.find_cluster(cluster_id)
            } else {
                None
            };
            Record {
                align,
                domain,
                cluster,
            }
        })
        .collect()
}

/// `TDomSeqAlignIndex::__TCdAlignRecord` of the site output: a hit whose model's sites are
/// printed, with the model of its hierarchy root (`pRootCdInfo`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct FeatureRecord {
    pub(crate) align: usize,
    pub(crate) domain: usize,
    /// The root of the model's hierarchy, `None` when the model has none (or it is not in the
    /// data).
    pub(crate) root: Option<usize>,
}

/// `TDomSeqAlignIndex::CreateRecordSets`, the feature records (`rFeatAligns`) of `mode`: for
/// each representative hit, the qualified curated hits of its region (`ExtractFeatAligns`),
/// then the structural-domain hits.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1918-2032
// (The `idxidx` members only order the JSON output.)
pub(crate) fn feature_records(
    index: &AlignIndex,
    aligns: &[DomAlign],
    data: &CddData,
    mode: DataMode,
) -> Vec<FeatureRecord> {
    let indices = match mode {
        DataMode::Rep => &index.concise,
        DataMode::Std => &index.standard,
        DataMode::Full => &index.sorted,
    };
    let mut records = Vec::new();
    for &align in indices {
        let a = &aligns[align];
        if a.rep && a.rep_class == 0 {
            let domain = data
                .find_domain(a.pssm_id)
                .expect("a classified hit has its model");
            extract_feat_aligns(align, domain, index, aligns, data, &mut records);
        }
    }
    for &align in &index.structural {
        let domain = data
            .find_domain(aligns[align].pssm_id)
            .expect("a classified hit has its model");
        records.push(FeatureRecord {
            align,
            domain,
            root: None,
        });
    }
    records
}

/// `TDomSeqAlignIndex::ExtractFeatAligns`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:2208-2272
// `rec` is one object for the whole call: a model without a hierarchy root (`m_uiHierarchyRoot`
// 0) keeps the `pRootCdInfo` of the record pushed before it (`root` below, set only when the
// model has a root).
fn extract_feat_aligns(
    rep_align: usize,
    rep_domain: usize,
    index: &AlignIndex,
    aligns: &[DomAlign],
    data: &CddData,
    result: &mut Vec<FeatureRecord>,
) {
    let rep_cluster = data.domains[rep_domain].cluster_id;
    let mut root: Option<usize> = None;
    for &candidate in &index.qualified_feat {
        // One model may have several hits; each hit is listed once.
        if result.iter().any(|record| record.align == candidate) {
            continue;
        }
        if aligns[candidate].region != aligns[rep_align].region {
            continue;
        }
        let domain = data
            .find_domain(aligns[candidate].pssm_id)
            .expect("a classified hit has its model");
        if rep_cluster > 0 && data.domains[domain].cluster_id != rep_cluster {
            continue;
        }
        let hierarchy_root = data.domains[domain].hierarchy_root;
        if hierarchy_root > 0 {
            root = if hierarchy_root == aligns[candidate].pssm_id {
                Some(domain)
            } else {
                data.find_domain(hierarchy_root)
            };
        }
        result.push(FeatureRecord {
            align: candidate,
            domain,
            root,
        });
    }
}

#[cfg(test)]
#[path = "domains_tests.rs"]
mod tests;
