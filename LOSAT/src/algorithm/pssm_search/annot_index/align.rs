//! The rows of a Seq-align as NCBI's object manager indexes them
//! (`CAnnotObject_Info::x_ProcessAlign`): one `HandleRangeMap` per row, with the warnings of the
//! sizes that do not agree.

use super::range::{
    interval_parts, point_parts, Handle, HandleRangeMap, LiveHandles, Range, Strand, UNKNOWN,
};
use super::view::V;
use super::{Indexer, Refusal, Stop};
use crate::ncbi_exception::{NcbiException, Severity, ThrowSite};

/// Std-seg's `dim` above which NCBI's rows depend on the machine's memory (`hrmaps.resize(dim)`
/// allocates `dim` maps; line RG D3-1 `om/SPEC.md` §10 Q2).
const STD_DIM_LIMIT: u64 = 1_000_000;

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:1548-1552
// ```c++
//     default:
//         NCBI_THROW_FMT(CSeqLocException, eUnsupported,
//                        "CSeq_loc_CI: unsupported location type: "<<
//                        loc.SelectionName(loc.Which()));
// ```
const LOC_CI_UNSUPPORTED: ThrowSite = ThrowSite {
    file: "src/objects/seq/../seqloc/Seq_loc.cpp",
    line: 1552,
    kind: "CSeqLocException::eUnsupported",
    function: Some("ncbi::objects::CSeq_loc_CI_Impl::x_ProcessLocation"),
};

/// The element `index` of a list member that the checks before keep in range.
fn at(list: V<'_>, index: usize) -> i64 {
    list.int_in(index)
}

// NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:563-846 (`x_ProcessAlign`, quoted per
// segment type below)
/// `x_ProcessAlign(rows, align)`: the rows of `align` added to `rows`, with NCBI's warnings
/// (`Warning: (1207.<n>)`).
pub(super) fn process_align(
    ix: &mut Indexer<'_>,
    rows: &mut Vec<HandleRangeMap>,
    align: V<'_>,
) -> Result<(), Stop> {
    let Some(segs) = align.get("segs") else {
        return Ok(());
    };
    let (variant, segs) = segs.choice();
    match variant {
        "dendiag" => {
            for diag in segs.items() {
                dense_diag(ix, rows, diag);
            }
            Ok(())
        }
        "denseg" => dense_seg(ix, rows, segs),
        "std" => {
            for seg in segs.items() {
                std_seg(&mut ix.handles, rows, seg)?;
            }
            Ok(())
        }
        "packed" => packed_seg(&mut ix.handles, rows, segs),
        "disc" => {
            for sub in segs.items() {
                process_align(ix, rows, sub)?;
            }
            Ok(())
        }
        "spliced" => {
            spliced_seg(&mut ix.handles, rows, segs);
            Ok(())
        }
        "sparse" => {
            sparse_seg(ix, rows, segs);
            Ok(())
        }
        other => unreachable!("Seq-align.segs has no variant {other}"),
    }
}

/// `hrmaps.resize(size)` that only grows (the callers' `if (hrmaps.size() < dim)`).
fn grow(rows: &mut Vec<HandleRangeMap>, size: usize) {
    if rows.len() < size {
        rows.resize_with(size, HandleRangeMap::default);
    }
}

/// `hrmaps[row]` of a row that the sizes made (`grow`).
fn row_map(rows: &mut [HandleRangeMap], row: usize) -> &mut HandleRangeMap {
    &mut rows[row]
}

// NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:571-604
// ```c++
//             ITERATE ( CSeq_align::C_Segs::TDendiag, it, dendiag ) {
//                 const CDense_diag& diag = **it;
//                 int dim = diag.GetDim();
//                 if (dim != (int)diag.GetIds().size()) {
//                     ERR_POST_X(1, Warning << "Invalid 'ids' size in dendiag");
//                     dim = min(dim, (int)diag.GetIds().size());
//                 }
//                 if (dim != (int)diag.GetStarts().size()) {
//                     ERR_POST_X(2, Warning << "Invalid 'starts' size in dendiag");
//                     dim = min(dim, (int)diag.GetStarts().size());
//                 }
//                 if (diag.IsSetStrands()
//                     && dim != (int)diag.GetStrands().size()) {
//                     ERR_POST_X(3, Warning << "Invalid 'strands' size in dendiag");
//                     dim = min(dim, (int)diag.GetStrands().size());
//                 }
//                 if ((int)hrmaps.size() < dim) {
//                     hrmaps.resize(dim);
//                 }
//                 TSeqPos len = (*it)->GetLen();
//                 for (int row = 0; row < dim; ++row) {
//                     const CSeq_id& id = *(*it)->GetIds()[row];
//                     TSeqPos from = (*it)->GetStarts()[row];
//                     TSeqPos to = from + len - 1;
//                     ENa_strand strand = eNa_strand_unknown;
//                     if ( (*it)->IsSetStrands() ) {
//                         strand = (*it)->GetStrands()[row];
//                     }
//                     hrmaps[row].SetMasterSeq(master);
//                     hrmaps[row].AddRange(id, from, to, strand);
//                 }
//             }
// ```
// `dim` is `int` (DEFAULT 2, Dense_diag_.hpp), compared with the sizes as `int`.
fn dense_diag(ix: &mut Indexer<'_>, rows: &mut Vec<HandleRangeMap>, diag: V<'_>) {
    let mut dim = diag.get("dim").map_or(2, V::int) as i32;
    let ids = diag.field("ids");
    let starts = diag.field("starts");
    let strands = diag.get("strands");
    if dim != ids.len() as i32 {
        ix.post(
            Severity::Warning,
            "1207.1",
            b"Invalid 'ids' size in dendiag",
        );
        dim = dim.min(ids.len() as i32);
    }
    if dim != starts.len() as i32 {
        ix.post(
            Severity::Warning,
            "1207.2",
            b"Invalid 'starts' size in dendiag",
        );
        dim = dim.min(starts.len() as i32);
    }
    if let Some(strands) = strands.filter(|strands| dim != strands.len() as i32) {
        ix.post(
            Severity::Warning,
            "1207.3",
            b"Invalid 'strands' size in dendiag",
        );
        dim = dim.min(strands.len() as i32);
    }
    let dim = dim.max(0) as usize;
    grow(rows, dim);
    let len = diag.field("len").int() as u32;
    let ids: Vec<V<'_>> = ids.items().collect();
    for (row, id) in ids.iter().enumerate().take(dim) {
        let from = at(starts, row) as u32;
        let to = from.wrapping_add(len).wrapping_sub(1);
        let strand = strands.map_or(UNKNOWN, |strands| at(strands, row));
        row_map(rows, row).add_range(id.seq_id(), Range::new(from, to), strand, &mut ix.handles);
    }
}

// NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:606-649
// ```c++
//             size_t dim    = size_t(denseg.GetDim());
//             size_t numseg = size_t(denseg.GetNumseg());
//             // claimed dimension may not be accurate :-/
//             if (numseg != denseg.GetLens().size()) {
//                 ERR_POST_X(4, Warning << "Invalid 'lens' size in denseg");
//                 numseg = min(numseg, denseg.GetLens().size());
//             }
//             if (dim != denseg.GetIds().size()) {
//                 ERR_POST_X(5, Warning << "Invalid 'ids' size in denseg");
//                 dim = min(dim, denseg.GetIds().size());
//             }
//             if (dim*numseg != denseg.GetStarts().size()) {
//                 ERR_POST_X(6, Warning << "Invalid 'starts' size in denseg");
//                 dim = min(dim*numseg, denseg.GetStarts().size()) / numseg;
//             }
//             if (denseg.IsSetStrands()
//                 && dim*numseg != denseg.GetStrands().size()) {
//                 ERR_POST_X(7, Warning << "Invalid 'strands' size in denseg");
//                 dim = min(dim*numseg, denseg.GetStrands().size()) / numseg;
//             }
//             if (hrmaps.size() < dim) {
//                 hrmaps.resize(dim);
//             }
//             for (size_t seg = 0;  seg < numseg;  seg++) {
//                 for (size_t row = 0;  row < dim;  row++) {
//                     if (denseg.GetStarts()[seg*dim + row] < 0 ) {
//                         continue;
//                     }
//                     const CSeq_id& id = *denseg.GetIds()[row];
//                     TSeqPos from = denseg.GetStarts()[seg*dim + row];
//                     TSeqPos to = from + denseg.GetLens()[seg] - 1;
//                     ENa_strand strand = eNa_strand_unknown;
//                     if ( denseg.IsSetStrands() ) {
//                         strand = denseg.GetStrands()[seg*dim + row];
//                     }
//                     hrmaps[row].SetMasterSeq(master);
//                     hrmaps[row].AddRange(id, from, to, strand);
//                 }
//             }
// ```
// `dim` and `numseg` are `int`s made `size_t` (a negative one is a huge size); after the first two
// checks both are below the sizes of `lens` and `ids`, so the products do not wrap. A `numseg` of
// 0 at a division is an integer division by zero (SIGFPE, exit status 136 after the warnings
// posted before it: probes `a07`, `a08`, `k11`): refused.
fn dense_seg(ix: &mut Indexer<'_>, rows: &mut Vec<HandleRangeMap>, seg: V<'_>) -> Result<(), Stop> {
    let mut dim = seg.get("dim").map_or(2, V::int) as i32 as i64 as u64;
    let mut numseg = seg.field("numseg").int() as i32 as i64 as u64;
    let ids = seg.field("ids");
    let starts = seg.field("starts");
    let lens = seg.field("lens");
    let strands = seg.get("strands");
    if numseg != lens.len() as u64 {
        ix.post(
            Severity::Warning,
            "1207.4",
            b"Invalid 'lens' size in denseg",
        );
        numseg = numseg.min(lens.len() as u64);
    }
    if dim != ids.len() as u64 {
        ix.post(Severity::Warning, "1207.5", b"Invalid 'ids' size in denseg");
        dim = dim.min(ids.len() as u64);
    }
    if dim * numseg != starts.len() as u64 {
        ix.post(
            Severity::Warning,
            "1207.6",
            b"Invalid 'starts' size in denseg",
        );
        if numseg == 0 {
            return Err(Stop::Refuse(Refusal::DivideByZero));
        }
        dim = (dim * numseg).min(starts.len() as u64) / numseg;
    }
    if let Some(strands) = strands.filter(|strands| dim * numseg != strands.len() as u64) {
        ix.post(
            Severity::Warning,
            "1207.7",
            b"Invalid 'strands' size in denseg",
        );
        if numseg == 0 {
            return Err(Stop::Refuse(Refusal::DivideByZero));
        }
        dim = (dim * numseg).min(strands.len() as u64) / numseg;
    }
    let (dim, numseg) = (dim as usize, numseg as usize);
    grow(rows, dim);
    let ids: Vec<V<'_>> = ids.items().collect();
    for seg_index in 0..numseg {
        let len = at(lens, seg_index) as u32;
        for (row, id) in ids.iter().enumerate().take(dim) {
            let start = at(starts, seg_index * dim + row);
            if start < 0 {
                continue;
            }
            let from = start as u32;
            let to = from.wrapping_add(len).wrapping_sub(1);
            let strand = strands.map_or(UNKNOWN, |strands| at(strands, seg_index * dim + row));
            row_map(rows, row).add_range(
                id.seq_id(),
                Range::new(from, to),
                strand,
                &mut ix.handles,
            );
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:651-678
// ```c++
//             ITERATE ( CSeq_align::C_Segs::TStd, it, std ) {
//                 size_t dim = size_t((*it)->GetDim());
//                 if (hrmaps.size() < dim) {
//                     hrmaps.resize((*it)->GetDim());
//                 }
//                 ITERATE ( CStd_seg::TLoc, it_loc, (*it)->GetLoc() ) {
//                     CSeq_loc_CI row_it(**it_loc);
//                     for (size_t row = 0; row_it; ++row_it, ++row) {
//                         if (row >= hrmaps.size()) {
//                             hrmaps.resize(row + 1);
//                         }
//                         const CSeq_id& id = row_it.GetSeq_id();
//                         TSeqPos from = row_it.GetRange().GetFrom();
//                         TSeqPos to = row_it.GetRange().GetTo();
//                         ENa_strand strand = row_it.GetStrand();
//                         hrmaps[row].SetMasterSeq(master);
//                         hrmaps[row].AddRange(id, from, to, strand);
//                     }
//                 }
//             }
// ```
// A negative `dim` is a `size_t` above `max_size()`: `resize` throws libstdc++'s
// `std::length_error("vector::_M_default_append")`, which `CATCH_ALL` reports (probe `a20`).
fn std_seg(
    handles: &mut LiveHandles,
    rows: &mut Vec<HandleRangeMap>,
    seg: V<'_>,
) -> Result<(), Stop> {
    let dim = seg.get("dim").map_or(2, V::int) as i32;
    if dim < 0 {
        return Err(Stop::std_exception("vector::_M_default_append"));
    }
    if rows.len() < dim as usize {
        if dim as u64 > STD_DIM_LIMIT {
            return Err(Stop::Refuse(Refusal::StdDim));
        }
        grow(rows, dim as usize);
    }
    for loc in seg.field("loc").items() {
        // The iterator's ranges are made when it is constructed (`x_ProcessLocation` of the whole
        // location), so a `feat` anywhere stops before any row.
        let mut ranges = Vec::new();
        let mut lookups = Vec::new();
        loc_ci_ranges(loc, &mut ranges, &mut lookups)?;
        // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:1420-1425 (the iterator
        // looks up the handle of every part when it is made, and `AddRange(id, ...)` looks it up
        // again; the iterator's handles go with it, after the location's rows)
        // ```c++
        // void CSeq_loc_CI_Impl::x_SetId(SSeq_loc_CI_RangeInfo& info,
        //                                const CSeq_id& id)
        // {
        //     info.m_Id = &id;
        //     info.m_IdHandle = CSeq_id_Handle::GetHandle(id);
        // }
        // ```
        for id in lookups {
            Handle::look_up_in_iterator(id, handles);
        }
        for (row, (id, range, strand)) in ranges.into_iter().enumerate() {
            grow(rows, row + 1);
            row_map(rows, row).add_range(id, range, strand, handles);
        }
        handles.end_iterator();
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:1407-1419 (`eEmpty_Skip`,
// `eOrder_Biological`: no reversal)
// ```c++
// CSeq_loc_CI_Impl::CSeq_loc_CI_Impl(const CSeq_loc&           loc,
//                                    CSeq_loc_CI::EEmptyFlag   empty_flag,
//                                    CSeq_loc_CI::ESeqLocOrder order)
//     ...
// {
//     x_ProcessLocation(loc);
//     if ( order == CSeq_loc_CI::eOrder_Positional  &&  loc.IsReverseStrand() ) {
//         reverse(m_Ranges.begin(), m_Ranges.end());
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:1439-1553 (`x_ProcessLocation`:
// `null`, `empty` (and not set) skipped, `whole` the whole range, `int`, `pnt`, `packed-int`,
// `packed-pnt` (one id and strand), `mix`, `equiv` (every part), `bond` (a, then b), and the
// exception above for any other)
/// The ranges of `CSeq_loc_CI(loc)`: Seq-id, range and strand (`unknown` when not set); and the
/// Seq-ids its construction looks up (`x_SetId`: one for each part, one for a `packed-pnt`).
fn loc_ci_ranges<'a>(
    loc: V<'a>,
    out: &mut Vec<(&'a crate::seqid::ParsedSeqId, Range, Strand)>,
    lookups: &mut Vec<&'a crate::seqid::ParsedSeqId>,
) -> Result<(), Stop> {
    let (variant, value) = loc.choice();
    let mut push = |part: (&'a crate::seqid::ParsedSeqId, Range, Strand)| {
        lookups.push(part.0);
        out.push(part);
    };
    match variant {
        "null" | "empty" => {}
        "whole" => push((value.seq_id(), Range::WHOLE, UNKNOWN)),
        "int" => push(interval_parts(value)),
        "pnt" => push(point_parts(value)),
        "packed-int" => {
            for part in value.items().map(interval_parts) {
                push(part);
            }
        }
        "packed-pnt" => {
            let id = value.field("id").seq_id();
            let strand = value.get("strand").map_or(UNKNOWN, V::int);
            lookups.push(id);
            for point in value.field("points").ints() {
                let point = point as u32;
                out.push((id, Range::new(point, point), strand));
            }
        }
        "mix" | "equiv" => {
            for part in value.items() {
                loc_ci_ranges(part, out, lookups)?;
            }
        }
        "bond" => {
            push(point_parts(value.field("a")));
            if let Some(b) = value.get("b") {
                push(point_parts(b));
            }
        }
        other => {
            return Err(Stop::exception(NcbiException::new(
                LOC_CI_UNSUPPORTED,
                format!("CSeq_loc_CI: unsupported location type: {other}"),
            )))
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:680-716
// ```c++
//             size_t dim    = size_t(packed.GetDim());
//             size_t numseg = size_t(packed.GetNumseg());
//             // claimed dimension may not be accurate :-/
//             if (dim * numseg > packed.GetStarts().size()) {
//                 dim = packed.GetStarts().size() / numseg;
//             }
//             if (dim * numseg > packed.GetPresent().size()) {
//                 dim = packed.GetPresent().size() / numseg;
//             }
//             if (dim > packed.GetLens().size()) {
//                 dim = packed.GetLens().size();
//             }
//             if (hrmaps.size() < dim) {
//                 hrmaps.resize(dim);
//             }
//             for (size_t seg = 0;  seg < numseg;  seg++) {
//                 for (size_t row = 0;  row < dim;  row++) {
//                     if ( packed.GetPresent()[seg*dim + row] ) {
//                         hrmaps[row].SetMasterSeq(master);
//                         const CSeq_id& id = *packed.GetIds()[row];
//                         TSeqPos from = packed.GetStarts()[seg*dim + row];
//                         TSeqPos to = from + packed.GetLens()[seg] - 1;
//                         ENa_strand strand = eNa_strand_unknown;
//                         if ( packed.IsSetStrands() ) {
//                             strand = packed.GetStrands()[seg*dim + row];
//                         }
//                         hrmaps[row].AddRange(id, from, to, strand);
//                     }
//                 }
//             }
// ```
// The products are `size_t` and wrap. `ids` (a row past its size) and `lens` (a segment past its
// size) are read without a check: NCBI reads past them (refused). A strand read past `strands`
// changes no output of the index (`annot_index::index_align`: an alignment's row only posts when
// all its ranges are empty, which no strand changes), so it is read as `unknown`.
fn packed_seg(
    handles: &mut LiveHandles,
    rows: &mut Vec<HandleRangeMap>,
    seg: V<'_>,
) -> Result<(), Stop> {
    let mut dim = seg.get("dim").map_or(2, V::int) as i32 as i64 as u64;
    let numseg = seg.field("numseg").int() as i32 as i64 as u64;
    let ids = seg.field("ids");
    let starts = seg.field("starts");
    let present = seg.field("present").octets();
    let lens = seg.field("lens");
    let strands = seg.get("strands");
    if dim.wrapping_mul(numseg) > starts.len() as u64 {
        dim = starts.len() as u64 / numseg;
    }
    if dim.wrapping_mul(numseg) > present.len() as u64 {
        dim = present.len() as u64 / numseg;
    }
    if dim > lens.len() as u64 {
        dim = lens.len() as u64;
    }
    let dim = dim as usize;
    grow(rows, dim);
    if dim == 0 {
        return Ok(());
    }
    // A `numseg` whose product with `dim` wraps reads `present` and `starts` past their ends, as
    // do `ids` and `lens` shorter than `dim` and `numseg`: refused when reached. The loop stops at
    // the first read past `present`, so it is bounded by its size.
    let ids: Vec<V<'_>> = ids.items().collect();
    for seg_index in 0..numseg {
        for row in 0..dim {
            let i = (seg_index as usize)
                .checked_mul(dim)
                .and_then(|base| base.checked_add(row));
            let Some(&flag) = i.and_then(|i| present.get(i)) else {
                return Err(Stop::Refuse(Refusal::PackedOutOfBounds));
            };
            if flag == 0 {
                continue;
            }
            let i = i.expect("checked above");
            let (Some(id), Some(len), Some(start)) = (
                ids.get(row),
                lens.int_at(seg_index as usize),
                starts.int_at(i),
            ) else {
                return Err(Stop::Refuse(Refusal::PackedOutOfBounds));
            };
            let from = start as u32;
            let to = from.wrapping_add(len as u32).wrapping_sub(1);
            let strand = strands
                .and_then(|strands| strands.int_at(i))
                .unwrap_or(UNKNOWN);
            row_map(rows, row).add_range(id.seq_id(), Range::new(from, to), strand, handles);
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:727-778
// ```c++
//             const CSeq_id* gen_id = spliced.IsSetGenomic_id() ?
//                 &spliced.GetGenomic_id() : 0;
//             const CSeq_id* prod_id = spliced.IsSetProduct_id() ?
//                 &spliced.GetProduct_id() : 0;
//             hrmaps.resize(2);
//             ITERATE ( CSpliced_seg::TExons, it, spliced.GetExons() ) {
//                 const CSpliced_exon& ex = **it;
//                 const CSeq_id* ex_gen_id = ex.IsSetGenomic_id() ?
//                     &ex.GetGenomic_id() : gen_id;
//                 if ( ex_gen_id ) {
//                     ...
//                     TSeqPos from = ex.GetGenomic_start();
//                     TSeqPos to = ex.GetGenomic_end();
//                     ENa_strand strand = eNa_strand_unknown;
//                     if ( ex.IsSetGenomic_strand() ) {
//                         strand = ex.GetGenomic_strand();
//                     }
//                     else if ( spliced.IsSetGenomic_strand() ) {
//                         strand = spliced.GetGenomic_strand();
//                     }
//                     hrmaps[1].SetMasterSeq(master);
//                     hrmaps[1].AddRange(id, from, to, strand);
//                 }
//                 const CSeq_id* ex_prod_id = ex.IsSetProduct_id() ?
//                     &ex.GetProduct_id() : prod_id;
//                 if ( ex_prod_id ) {
//                     const CSeq_id& id = *ex_prod_id;
//                     TSeqPos from =
//                         (ex.GetProduct_start().IsNucpos() ?
//                          ex.GetProduct_start().GetNucpos() :
//                          ex.GetProduct_start().GetProtpos().GetAmin());
//                     TSeqPos to = ... (the same of GetProduct_end())
//                     ENa_strand strand = ... (product-strand of the exon, else of the spliced seg)
//                     hrmaps[0].SetMasterSeq(master);
//                     hrmaps[0].AddRange(id, from, to, strand);
//                 }
//             }
// ```
// `resize(2)` also cuts rows that a Seq-align before it in the same `disc` made.
fn spliced_seg(handles: &mut LiveHandles, rows: &mut Vec<HandleRangeMap>, seg: V<'_>) {
    rows.truncate(2);
    grow(rows, 2);
    let genomic_id = seg.get("genomic-id");
    let product_id = seg.get("product-id");
    let product_pos = |pos: V<'_>| -> u32 {
        let (variant, value) = pos.choice();
        match variant {
            "nucpos" => value.int() as u32,
            _ => value.field("amin").int() as u32,
        }
    };
    for exon in seg.field("exons").items() {
        if let Some(id) = exon.get("genomic-id").or(genomic_id) {
            let from = exon.field("genomic-start").int() as u32;
            let to = exon.field("genomic-end").int() as u32;
            let strand = exon
                .get("genomic-strand")
                .or(seg.get("genomic-strand"))
                .map_or(UNKNOWN, V::int);
            row_map(rows, 1).add_range(id.seq_id(), Range::new(from, to), strand, handles);
        }
        if let Some(id) = exon.get("product-id").or(product_id) {
            let from = product_pos(exon.field("product-start"));
            let to = product_pos(exon.field("product-end"));
            let strand = exon
                .get("product-strand")
                .or(seg.get("product-strand"))
                .map_or(UNKNOWN, V::int);
            row_map(rows, 0).add_range(id.seq_id(), Range::new(from, to), strand, handles);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:780-836
// ```c++
//             // consensus sequence row + one row for each sub-alignment
//             size_t dim = sparse.GetRows().size() + 1;
//             if (hrmaps.size() < dim) {
//                 hrmaps.resize(dim);
//             }
//             size_t row = 0;
//             hrmaps[0].SetMasterSeq(master);
//             ITERATE ( CSparse_seg::TRows, it, sparse.GetRows() ) {
//                 const CSparse_align& aln_row = **it;
//                 size_t numseg = aln_row.GetNumseg();
//                 if (numseg != aln_row.GetFirst_starts().size()) {
//                     ERR_POST_X(9, Warning <<
//                         "Invalid size of 'first-starts' in sparse-align");
//                     numseg = min(numseg, aln_row.GetFirst_starts().size());
//                 }
//                 if (numseg != aln_row.GetSecond_starts().size()) {
//                     ERR_POST_X(10, Warning <<
//                         "Invalid size of 'second-starts' in sparse-align");
//                     numseg = min(numseg, aln_row.GetSecond_starts().size());
//                 }
//                 if (numseg != aln_row.GetLens().size()) {
//                     ERR_POST_X(11, Warning <<
//                         "Invalid size of 'lens' in sparse-align");
//                     numseg = min(numseg, aln_row.GetLens().size());
//                 }
//                 if (aln_row.IsSetSecond_strands()  &&
//                     numseg != aln_row.GetSecond_strands().size()) {
//                     ERR_POST_X(12, Warning <<
//                         "Invalid size of 'second-strands' in sparse-align");
//                     numseg = min(numseg, aln_row.GetSecond_strands().size());
//                 }
//
//                 hrmaps[row+1].SetMasterSeq(master);
//                 for (size_t seg = 0; seg < numseg; ++seg) {
//                     TSeqPos len = aln_row.GetLens()[seg];
//                     {
//                         const CSeq_id& id = aln_row.GetFirst_id();
//                         TSeqPos from = aln_row.GetFirst_starts()[seg];
//                         TSeqPos to = from + len - 1;
//                         // consensus sequence goes to the first row
//                         hrmaps[0].AddRange(id, from, to);
//                     }
//                     {
//                         const CSeq_id& id = aln_row.GetSecond_id();
//                         TSeqPos from = aln_row.GetSecond_starts()[seg];
//                         TSeqPos to = from + len - 1;
//                         ENa_strand strand = eNa_strand_unknown;
//                         if ( aln_row.IsSetSecond_strands() ) {
//                             strand = aln_row.GetSecond_strands()[row];
//                         }
//                         hrmaps[row+1].AddRange(id, from, to, strand);
//                     }
//                 }
//                 row++;
//             }
// ```
// The second strand is read at the row's index, not the segment's (past `second-strands` for a
// row beyond its size): no strand changes an output of the index (`packed_seg`), so it is read as
// `unknown` there.
fn sparse_seg(ix: &mut Indexer<'_>, rows: &mut Vec<HandleRangeMap>, seg: V<'_>) {
    let sparse_rows: Vec<V<'_>> = seg.field("rows").items().collect();
    grow(rows, sparse_rows.len() + 1);
    for (row, aln) in sparse_rows.iter().enumerate() {
        let mut numseg = aln.field("numseg").int() as i32 as i64 as u64;
        let first_starts = aln.field("first-starts");
        let second_starts = aln.field("second-starts");
        let lens = aln.field("lens");
        let second_strands = aln.get("second-strands");
        if numseg != first_starts.len() as u64 {
            ix.post(
                Severity::Warning,
                "1207.9",
                b"Invalid size of 'first-starts' in sparse-align",
            );
            numseg = numseg.min(first_starts.len() as u64);
        }
        if numseg != second_starts.len() as u64 {
            ix.post(
                Severity::Warning,
                "1207.10",
                b"Invalid size of 'second-starts' in sparse-align",
            );
            numseg = numseg.min(second_starts.len() as u64);
        }
        if numseg != lens.len() as u64 {
            ix.post(
                Severity::Warning,
                "1207.11",
                b"Invalid size of 'lens' in sparse-align",
            );
            numseg = numseg.min(lens.len() as u64);
        }
        if let Some(strands) = second_strands.filter(|s| numseg != s.len() as u64) {
            ix.post(
                Severity::Warning,
                "1207.12",
                b"Invalid size of 'second-strands' in sparse-align",
            );
            numseg = numseg.min(strands.len() as u64);
        }
        let first_id = aln.field("first-id").seq_id();
        let second_id = aln.field("second-id").seq_id();
        for seg_index in 0..numseg as usize {
            let len = at(lens, seg_index) as u32;
            let from = at(first_starts, seg_index) as u32;
            row_map(rows, 0).add_range(
                first_id,
                Range::new(from, from.wrapping_add(len).wrapping_sub(1)),
                UNKNOWN,
                &mut ix.handles,
            );
            let from = at(second_starts, seg_index) as u32;
            let strand = second_strands
                .and_then(|strands| strands.int_at(row))
                .unwrap_or(UNKNOWN);
            row_map(rows, row + 1).add_range(
                second_id,
                Range::new(from, from.wrapping_add(len).wrapping_sub(1)),
                strand,
                &mut ix.handles,
            );
        }
    }
}
