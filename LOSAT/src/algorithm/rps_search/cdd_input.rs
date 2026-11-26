//! The conserved-domain input of DELTA-BLAST's PSSM (`CCddInputData`): the domain hits of the
//! RPS search are filtered by their E-value, the overlapping hits of one domain are trimmed, and
//! the alignment of the hits' columns with the query is made, which `core::psi::cd_stages` turns
//! into the PSSM.
//!
//! NCBI keeps the ranges of a hit as `CRange<int>` made of a start and an end that the range
//! reads as closed (`GetTo` is the last position, `GetToOpen` one past it), while `CCddInputData`
//! builds them as start and start plus length: a segment of `n` columns at `f` is the range
//! `f..f+n`, `n + 1` positions long. The intersection, the emptiness and the comparisons of the
//! trimming read that range, the columns of the alignment its `n`. `Range` holds `from` and the
//! closed `to` so that each test below is NCBI's.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:96-330
//! NCBI reference (598d8ae6): c++/include/algo/blast/api/cdd_pssm_input.hpp:195-400
//! NCBI reference (598d8ae6): c++/include/util/range.hpp:105-381

use std::cmp::Ordering;

use crate::core::psi::{CdCell, CdColumn, PsiCdMsa};
use crate::seqdb::rps::{CdProfileFiles, RpsOpenError};
use crate::seqid::model::{ObjectId, ParsedSeqId};

/// A `CRange<int>` as `CCddInputData` uses it: `from` and the closed end `to` (`GetTo`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Range {
    pub(crate) from: i32,
    pub(crate) to: i32,
}

impl Range {
    /// `COpenRange::Empty`: the open end `to + 1` is not past `from`.
    fn is_empty(self) -> bool {
        self.to < self.from
    }

    /// `COpenRange::IntersectionWith`: the later start and the earlier end.
    fn intersection(self, other: Range) -> Range {
        Range {
            from: self.from.max(other.from),
            to: self.to.min(other.to),
        }
    }
}

/// One aligned segment of a hit (`CHitSegment`): a run of the query and of the domain's columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Segment {
    pub(crate) query: Range,
    pub(crate) subject: Range,
}

impl Segment {
    /// `CHitSegment::AdjustRanges`.
    fn adjust(&mut self, d_from: i32, d_to: i32) {
        self.query.from += d_from;
        self.query.to += d_to;
        self.subject.from += d_from;
        self.subject.to += d_to;
    }

    /// `CHitSegment::IsEmpty`.
    fn is_empty(&self) -> bool {
        self.query.from > self.query.to || self.subject.from > self.subject.to
    }
}

/// A hit of the domain search (`CHit`): the domain, its E-value and the aligned segments.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Hit {
    /// The Seq-align id of the domain.
    pub(crate) id: ParsedSeqId,
    /// The domain's OID in the domain database, whose files hold its columns.
    pub(crate) oid: u32,
    pub(crate) e_value: f64,
    pub(crate) segments: Vec<Segment>,
}

/// What LOSAT does not run as NCBI does.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum CddRefusal {
    /// A domain whose Seq-align id is neither a local nor a general id: NCBI orders them by the
    /// rules of `CSeq_id::CompareOrdered` that LOSAT does not port.
    IdKind(String),
    /// A hit segment that reaches past the query or has negative length: NCBI writes outside its
    /// alignment.
    Range(String),
}

impl CddRefusal {
    /// The text after `... is not supported by LOSAT's DELTABLAST`'s subject.
    pub(crate) fn message(&self) -> String {
        match self {
            CddRefusal::IdKind(id) => format!(
                "a domain database whose Seq-align ids are not local or general ids (the hit {id}) is not supported by LOSAT's DELTABLAST"
            ),
            CddRefusal::Range(id) => format!(
                "a query whose domain hit {id} has a segment that is outside the query or has a negative length (NCBI BLAST+ writes outside its alignment) is not supported by LOSAT's DELTABLAST"
            ),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:335-376
/// The hit of a gapped alignment: its segments of aligned letters (the gaps of either sequence
/// are dropped), each `len` columns from `query_start` and `subject_start` (0-based).
pub(crate) fn hit_of_segments(
    id: ParsedSeqId,
    oid: u32,
    e_value: f64,
    aligned: impl IntoIterator<Item = (i32, i32, i32)>,
) -> Hit {
    Hit {
        id,
        oid,
        e_value,
        segments: aligned
            .into_iter()
            .map(|(query, subject, len)| Segment {
                query: Range {
                    from: query,
                    to: query + len,
                },
                subject: Range {
                    from: subject,
                    to: subject + len,
                },
            })
            .collect(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:151-165
/// `x_ProcessAlignments`: the hits whose E-value is below the domains' inclusion threshold (the
/// lowest accepted is -1).
pub(crate) fn included_hits(hits: Vec<Hit>, inclusion_ethresh: f64) -> Vec<Hit> {
    hits.into_iter()
        .filter(|hit| hit.e_value >= -1.0 && hit.e_value < inclusion_ethresh)
        .collect()
}

impl Hit {
    /// `CHit::IsEmpty`.
    fn is_empty(&self) -> bool {
        self.segments.iter().all(Segment::is_empty)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:442-538
    /// The intersection of `CHit` with a list of ranges: the segments cut to the sorted `ranges`,
    /// by the subject's positions or the query's.
    fn intersect_with_ranges(&mut self, ranges: &[Range], by_subject: bool) {
        let mut slots: Vec<Option<Segment>> = std::mem::take(&mut self.segments)
            .into_iter()
            .map(Some)
            .collect();
        let mut new_segments: Vec<Segment> = Vec::new();
        let mut range = 0;
        let mut at = 0;
        while at < slots.len() && range < ranges.len() {
            let segment = slots[at].expect("segments before `at` are looked at once");
            let seg_range = if by_subject {
                segment.subject
            } else {
                segment.query
            };
            // The ranges strictly below the segment.
            while range < ranges.len() && ranges[range].to < seg_range.from {
                range += 1;
            }
            if range == ranges.len() {
                break;
            }
            let intersection = seg_range.intersection(ranges[range]);
            if intersection == seg_range {
                at += 1;
                continue;
            }
            if intersection.is_empty() {
                slots[at] = None;
                at += 1;
                continue;
            }
            while range < ranges.len() && ranges[range].from < seg_range.to {
                let d_from = seg_range.from.max(ranges[range].from) - seg_range.from;
                let d_to = seg_range.to.min(ranges[range].to) - seg_range.to;
                let mut piece = segment;
                piece.adjust(d_from, d_to);
                new_segments.push(piece);
                range += 1;
            }
            slots[at] = None;
            at += 1;
        }
        // The segments behind the last range have an empty intersection.
        for slot in &mut slots[at..] {
            *slot = None;
        }
        new_segments.extend(slots.into_iter().flatten());
        // `compare_hitseg_range`: by the subject's start (the starts are distinct).
        new_segments.sort_by_key(|segment| segment.subject.from);
        self.segments = new_segments;
    }

    /// The intersection of `CHit` with another hit: this hit's segments cut to the other hit's ranges.
    fn intersect_with(&mut self, other: &Hit, by_subject: bool) {
        let mut ranges: Vec<Range> = other
            .segments
            .iter()
            .map(|segment| {
                if by_subject {
                    segment.subject
                } else {
                    segment.query
                }
            })
            .collect();
        // `compare_range`
        ranges.sort_by_key(|range| (range.from, range.to));
        self.intersect_with_ranges(&ranges, by_subject);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:569-755
    /// `CHit::Subtract`: takes the query range that `hit` covers out of this hit's segments.
    fn subtract(&mut self, hit: &Hit) {
        if self.is_empty() || hit.is_empty() {
            return;
        }
        let (Some(first), Some(last)) = (hit.segments.first(), hit.segments.last()) else {
            return;
        };
        let from = first.query.from;
        let to = last.query.to;
        let segments = &self.segments;
        // No overlap between the hits.
        if segments[0].query.from >= to || segments[segments.len() - 1].query.to <= from {
            return;
        }
        let mut new_segments: Vec<Segment> = Vec::with_capacity(segments.len());
        let mut at = 0;
        // The segments that end before the subtracted hit starts stay.
        while at < segments.len() && segments[at].query.to <= from {
            new_segments.push(segments[at]);
            at += 1;
        }
        // Nothing of the hit overlaps the subtracted one: no change.
        if at == segments.len() || segments[at].query.from > to {
            return;
        }
        if segments[at].query.to > to {
            // The segment covers the whole subtracted hit: what is left and right of it.
            let mut right = segments[at];
            if segments[at].query.from < from {
                let mut left = segments[at];
                left.adjust(0, from - segments[at].query.to);
                new_segments.push(left);
            }
            right.adjust(to - right.query.from, 0);
            new_segments.push(right);
            at += 1;
            new_segments.extend_from_slice(&segments[at..]);
        } else {
            // A segment inside the subtracted hit is deleted, one that starts before it ends there.
            if segments[at].query.from < from {
                let mut kept = segments[at];
                kept.adjust(0, from - segments[at].query.to);
                new_segments.push(kept);
            }
            at += 1;
            while at < segments.len() && segments[at].query.to <= to {
                at += 1;
            }
            if at < segments.len() {
                if segments[at].query.from < to {
                    let mut kept = segments[at];
                    kept.adjust(to - segments[at].query.from, 0);
                    new_segments.push(kept);
                }
                at += 1;
                new_segments.extend_from_slice(&segments[at..]);
            }
        }
        self.segments = new_segments;
    }
}

/// The kind of a Seq-id in `CSeq_id::Which` order, for the ids that LOSAT orders.
fn id_choice(id: &ParsedSeqId) -> Option<u8> {
    match id {
        ParsedSeqId::Local(_) => Some(1),
        ParsedSeqId::General(_) => Some(11),
        _ => None,
    }
}

/// `CObject_id::GetIdType` and `Compare`: a string that reads as an integer without a leading
/// zero or plus sign is the integer.
enum ObjectKey<'a> {
    Id(i64),
    Str(&'a [u8]),
}

fn object_key(tag: &ObjectId) -> ObjectKey<'_> {
    match tag {
        ObjectId::Id(id) => ObjectKey::Id(i64::from(*id)),
        ObjectId::Str(text) => {
            let digits = text.strip_prefix(b"-").unwrap_or(text);
            let canonical = !digits.is_empty()
                && digits.iter().all(u8::is_ascii_digit)
                && (digits.len() == 1 || digits[0] != b'0')
                && !(text.len() == 2 && text[0] == b'-' && digits == b"0");
            match std::str::from_utf8(text)
                .ok()
                .filter(|_| canonical)
                .and_then(|text| text.parse::<i64>().ok())
            {
                Some(value) => ObjectKey::Id(value),
                None => ObjectKey::Str(text),
            }
        }
    }
}

/// `PNocase().Compare`: case-insensitive comparison of the bytes, the shorter first on a common
/// start.
fn compare_nocase(a: &[u8], b: &[u8]) -> Ordering {
    let lower = |bytes: &[u8]| bytes.iter().map(u8::to_ascii_lowercase).collect::<Vec<_>>();
    lower(a).cmp(&lower(b))
}

fn compare_object(a: &ObjectId, b: &ObjectId) -> Ordering {
    match (object_key(a), object_key(b)) {
        (ObjectKey::Id(a), ObjectKey::Id(b)) => a.cmp(&b),
        (ObjectKey::Id(_), ObjectKey::Str(_)) => Ordering::Less,
        (ObjectKey::Str(_), ObjectKey::Id(_)) => Ordering::Greater,
        (ObjectKey::Str(a), ObjectKey::Str(b)) => compare_nocase(a, b),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:486-522
// NCBI reference (598d8ae6): c++/src/objects/general/Dbtag.cpp:176-183
// NCBI reference (598d8ae6): c++/src/objects/general/Object_id.cpp:145-161
/// `CSeq_id::CompareOrdered` of two local or general ids (`id_choice`).
fn compare_ordered(a: &ParsedSeqId, b: &ParsedSeqId) -> Ordering {
    match (a, b) {
        (ParsedSeqId::Local(a), ParsedSeqId::Local(b)) => compare_object(a, b),
        (ParsedSeqId::General(a), ParsedSeqId::General(b)) => {
            compare_nocase(&a.db, &b.db).then_with(|| compare_object(&a.tag, &b.tag))
        }
        _ => id_choice(a).cmp(&id_choice(b)),
    }
}

/// The text of an id for LOSAT's messages.
fn id_label(id: &ParsedSeqId) -> String {
    let mut out = Vec::new();
    crate::seqid::label::fasta(id, &mut out);
    String::from_utf8_lossy(&out).into_owned()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:168-219
/// `x_RemoveMultipleCdHits`: the hits sorted by domain (`CSeq_id::CompareOrdered`) and E-value,
/// and each trimmed by the kept hits of its domain so that no query range of a domain is counted
/// twice; a hit left empty is dropped.
pub(crate) fn remove_multiple_cd_hits(mut hits: Vec<Hit>) -> Result<Vec<Hit>, CddRefusal> {
    if hits.len() < 2 {
        return Ok(hits);
    }
    for hit in &hits {
        if id_choice(&hit.id).is_none() {
            return Err(CddRefusal::IdKind(id_label(&hit.id)));
        }
    }
    // `compare_hits_by_seqid_eval`: the domain's id, then the E-value. Rust's stable sort keeps hits
    // of one domain with one E-value in the order of the domain search's results; NCBI's sort
    // leaves them in the order its C++ library's sort gives, which can trim a hit of the domain
    // differently (PD-LOSAT-NCBI-DEFECTS, approved difference of version 1.116; Owner, 2026-11-24).
    hits.sort_by(|a, b| {
        compare_ordered(&a.id, &b.id)
            .then_with(|| a.e_value.partial_cmp(&b.e_value).unwrap_or(Ordering::Equal))
    });
    let mut kept: Vec<Hit> = Vec::with_capacity(hits.len());
    let mut hits = hits.into_iter();
    kept.extend(hits.next());
    for mut hit in hits {
        // Each kept hit of the same domain, the last first.
        let mut emptied = false;
        for earlier in kept.iter().rev() {
            if compare_ordered(&hit.id, &earlier.id) != Ordering::Equal {
                break;
            }
            // The kept hit cut to this hit's subject ranges, with its own query ranges.
            let mut intersection = earlier.clone();
            intersection.intersect_with(&hit, true);
            hit.subtract(&intersection);
            if hit.is_empty() {
                emptied = true;
                break;
            }
        }
        if !emptied {
            kept.push(hit);
        }
    }
    Ok(kept)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:222-272
/// `x_FillHitsData` and `x_CreateMsa`: the alignment of the kept hits' columns with the query,
/// one row for each hit with the columns of its segments in place; the columns' data come from the
/// domain database's `.wcounts` and `.obsr` files.
pub(crate) fn build_msa(
    query: &[u8],
    hits: &[Hit],
    files: &CdProfileFiles,
) -> Result<PsiCdMsa, anyhow::Error> {
    let mut columns: Vec<CdColumn> = Vec::new();
    let mut rows: Vec<Vec<CdCell>> = Vec::with_capacity(hits.len());
    for hit in hits {
        let mut row = vec![CdCell::UNALIGNED; query.len()];
        // The profile's columns that the hit's segments read, at most up to the last one's end
        // (`CdHitData`: a few reads for the hit).
        let end = hit
            .segments
            .iter()
            .map(|segment| {
                i64::from(segment.subject.from) + i64::from(segment.query.to - segment.query.from)
            })
            .max()
            .map_or(0, |end| usize::try_from(end).unwrap_or(0));
        let mut data = files.hit(hit.oid as usize, end);
        for segment in &hit.segments {
            let len = segment.query.to - segment.query.from;
            if len < 0
                || segment.query.from < 0
                || segment.subject.from < 0
                || segment.query.from as usize + len as usize > query.len()
            {
                anyhow::bail!(CddRefusal::Range(id_label(&hit.id)).message());
            }
            let from = segment.subject.from as usize;
            let data = data
                .columns(from, from + len as usize)
                .map_err(open_error)?;
            for (offset, (wfreqs, iobsr)) in data.wfreqs.into_iter().zip(data.iobsr).enumerate() {
                row[segment.query.from as usize + offset] = CdCell::aligned(columns.len());
                columns.push(CdColumn { wfreqs, iobsr });
            }
        }
        rows.push(row);
    }
    Ok(PsiCdMsa {
        query: query.to_vec(),
        rows,
        columns,
    })
}

/// NCBI's error of the RPS files as the engine's (exit status 3), LOSAT's refusal as its own.
pub(crate) fn open_error(error: RpsOpenError) -> anyhow::Error {
    match error {
        RpsOpenError::Ncbi(text) => crate::blastinput::app::engine_error(&text),
        RpsOpenError::Refused(text) => anyhow::anyhow!(text),
    }
}

#[cfg(test)]
#[path = "cdd_input_tests.rs"]
mod tests;
