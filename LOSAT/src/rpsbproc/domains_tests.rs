//! The pieces of the domain classification: `CleanAlignment`, `CalcMasterGaps`, `CountSrc`, the
//! `Std-seg`s of a translated query (`ParseAlignSegs`), its reading frames and the mapping of a
//! site to the frame's translation and to the plus strand (`GetTranslatedPosMap`,
//! `Pr2NaConvert`).

use super::*;
use crate::rpsbproc::archive::{Align, DenseSeg, Segs, StdLoc, StdPair};
use crate::rpsbproc::segset::Seg;

fn dense_seg(starts: &[i32], lens: &[u32]) -> DenseSeg {
    DenseSeg {
        ids: Vec::new(),
        starts: starts.to_vec(),
        lens: lens.to_vec(),
    }
}

fn align() -> Align {
    Align {
        segs: Segs::Other("dendiag"),
        score: Some(10),
        e_value: Some(1e-5),
        bit_score: Some(40.5),
        num_ident: None,
        score_of_other_type: None,
    }
}

// impl/basedata_impl.hpp:2-19 and biodata_blast.cpp:23-43 (quoted at `clean_alignment` and
// `parse_align`): the segments with a gap in either row are dropped; the range runs from the
// first kept segment to the end of the last one.
#[test]
fn parse_align_keeps_the_aligned_segments() {
    let seg = dense_seg(&[4, 0, -1, 70, 74, 80, 90, -1], &[70, 5, 10, 3]);
    let dom = parse_align(&seg, &align(), 42);
    assert_eq!(dom.m_starts, vec![4, 74]);
    assert_eq!(dom.s_starts, vec![0, 80]);
    assert_eq!(dom.lens, vec![70, 10]);
    assert_eq!((dom.from, dom.to), (4, 83));
    assert_eq!((dom.e_value, dom.bit_score, dom.pssm_id), (1e-5, 40.5, 42));
    // No aligned segment: NCBI's range stays `k_INT_MAX` .. 0.
    let none = parse_align(&dense_seg(&[-1, 0], &[5]), &align(), 42);
    assert_eq!((none.from, none.to), (i64::from(i32::MAX), 0));
    assert!(none.s_starts.is_empty());
}

// basedata.cpp:1324-1339 (quoted at `master_gaps`): query gaps of at least 35 residues.
#[test]
fn master_gaps_of_35_residues_and_more() {
    let seg = dense_seg(&[0, 0, 44, 10, 100, 20, 134, 30], &[10, 10, 10, 10]);
    let dom = parse_align(&seg, &align(), 1);
    let gaps = master_gaps(&dom);
    // 10..43 (34 residues) is not a gap; 54..99 (46) and 110..133 (24): only the first.
    assert_eq!(gaps.total_residues(), 46);
}

// basedata.cpp:154-166 (quoted at `count_source`): one hit per source database and region;
// the PRK family (PRK, CHL, MTH, PHA, PLN, PTZ) counts as one source; unknown sources always
// count.
#[test]
fn one_hit_per_source_database() {
    let mut counts = [0; 7];
    assert!(count_source(&mut counts, b"cd00001"));
    assert!(!count_source(&mut counts, b"cd00002"));
    assert!(count_source(&mut counts, b"pfam00001"));
    assert!(count_source(&mut counts, b"PRK00001"));
    assert!(!count_source(&mut counts, b"PTZ00001"));
    assert!(count_source(&mut counts, b"NF000001"));
    assert!(count_source(&mut counts, b"NF000002"));
}

fn int(from: u32, to: u32, strand: Option<i32>) -> StdLoc {
    StdLoc::Int { from, to, strand }
}

fn std_seg(query: StdLoc, hit: StdLoc) -> StdPair {
    StdPair { query, hit }
}

// biodata_blast.cpp:44-209 (quoted at `parse_std_align`): a plus-strand hit of rpstblastn
// (archive c01, Query_1: 387-2582 with 0-731, a query insertion 2583-2585 with an `empty`
// subject, 2586-2795 with 732-801): pseudo-protein starts `f / 3`, the frame of `from`.
#[test]
fn a_plus_strand_hit_of_a_translated_query() {
    let segs = [
        std_seg(int(387, 2582, Some(1)), int(0, 731, Some(0))),
        std_seg(int(2583, 2585, Some(1)), StdLoc::Empty),
        std_seg(int(2586, 2795, Some(1)), int(732, 801, Some(0))),
    ];
    let dom = parse_std_align(&segs, &align(), 7).unwrap();
    assert!(dom.pr2na);
    assert_eq!((dom.from, dom.to), (387, 2795));
    assert_eq!(dom.m_starts, vec![129, 862]);
    assert_eq!(dom.s_starts, vec![0, 732]);
    assert_eq!(dom.lens, vec![732, 70]);
    assert_eq!((dom.plus_rf_idx(), dom.aligned_len()), (0, 0));
    assert_eq!(rf_index(&dom, 5000), 0);
}

// The minus strand: `alignedLen = (m_iTo + 1) / 3`, starts `alignedLen - (f / 3 + len)`, and
// the frame index from the query's length (`PlusRFIdx2MinusRFIdx`).
#[test]
fn a_minus_strand_hit_of_a_translated_query() {
    let segs = [
        std_seg(int(301, 599, Some(2)), int(0, 99, None)),
        std_seg(int(100, 297, Some(2)), int(101, 166, None)),
    ];
    let dom = parse_std_align(&segs, &align(), 7).unwrap();
    assert_eq!((dom.from, dom.to), (100, 599));
    assert_eq!(dom.aligned_len(), 200);
    assert_eq!(dom.plus_rf_idx(), 1);
    // 200 - (301 / 3 + 100) = 0; 200 - (100 / 3 + 66) = 101.
    assert_eq!(dom.m_starts, vec![0, 101]);
    // (5000 - 1) % 3 + 3 = 4; (4999 - 1) % 3 + 3 = 3.
    assert_eq!(rf_index(&dom, 5000), 4);
    assert_eq!(rf_index(&dom, 4999), 3);
    let frames = sort_reading_frames(&[dom], 5000).unwrap();
    assert_eq!(frames[4], vec![0]);
    assert_eq!(
        sort_reading_frames(&[], 0),
        Err("Cannot sort reading frames - query sequence length unknown")
    );
}

// Without a strand the hit is an alignment of two proteins; a subject location that is not
// `int` or `pnt` leaves its query start alone (plus strand: kept as NCBI does), which the
// minus strand reads past (refused); `null` and `empty` segments are skipped; a `pnt` is a
// position of length 1. A subject interval whose length `t - f + 1` overflows `int` (undefined)
// has the length's value modulo 2^32, the `unsigned int` of `m_vecLens` (PD-LOSAT-NCBI-DEFECTS
// exception 11; PS16-Q1 oracle `ntslen` and its nearby input `ntslen1`, one position shorter).
#[test]
fn std_seg_forms() {
    let none =
        parse_std_align(&[std_seg(int(30, 59, None), int(0, 9, None))], &align(), 1).unwrap();
    assert!(!none.pr2na);
    assert_eq!(none.m_starts, vec![30]);
    let plus = [
        std_seg(int(30, 59, Some(1)), StdLoc::Whole),
        std_seg(StdLoc::Null, int(0, 9, None)),
        std_seg(
            StdLoc::Pnt {
                point: 90,
                strand: Some(3),
            },
            StdLoc::Pnt {
                point: 12,
                strand: None,
            },
        ),
    ];
    let dom = parse_std_align(&plus, &align(), 1).unwrap();
    assert_eq!(dom.m_starts, vec![10, 30]);
    assert_eq!(dom.s_starts, vec![12]);
    assert_eq!(dom.lens, vec![1]);
    let minus = [std_seg(int(30, 59, Some(2)), StdLoc::Whole)];
    assert_eq!(
        parse_std_align(&minus, &align(), 1).unwrap_err(),
        "a minus-strand Std-seg with a subject location other than int or pnt (NCBI reads past its segment lengths) is not supported by LOSAT's rpsbproc"
    );
    let length = |from: u32, to: u32| {
        let dom = parse_std_align(
            &[std_seg(int(0, 2, Some(1)), int(from, to, None))],
            &align(),
            1,
        )
        .unwrap();
        (dom.s_starts.get(0), dom.lens[0])
    };
    // 2^31 - 1 - (-2^31) + 1 = 2^32; 2^31 - 1 - (-5) + 1 = 2^31 + 5.
    assert_eq!(length(0x8000_0000, 0x7fff_ffff), (-2_147_483_648, 0));
    assert_eq!(length(0xffff_fffb, 0x7fff_ffff), (-5, 2_147_483_653));
    assert_eq!(length(0xffff_fffb, 0x7fff_fffe), (-5, 2_147_483_652));
    assert_eq!(length(0, 0x7fff_fffe), (0, 0x7fff_ffff));
}

// biodata_blast.cpp:63 (quoted at `parse_std_align`): `seqlocs[0]->IsNull() ||
// seqlocs[1]->IsNull()` reads the second location only when the first is not `null`. A
// segment whose only location is `null` is skipped (PS16c review oracle `nullone`: NCBI's
// output is that of the archive without the segment); a segment without locations, or whose
// only location is not `null`, is read past its end: refused.
#[test]
fn locations_are_read_in_ncbis_order() {
    let hit = std_seg(int(30, 59, Some(1)), int(0, 9, None));
    let alone = parse_std_align(&[hit], &align(), 1).unwrap();
    let with_null =
        parse_std_align(&[hit, std_seg(StdLoc::Null, StdLoc::Missing)], &align(), 1).unwrap();
    assert_eq!(
        (
            with_null.from,
            with_null.to,
            with_null.m_starts,
            with_null.lens
        ),
        (alone.from, alone.to, alone.m_starts, alone.lens)
    );
    let refusal = |segs: &[StdPair]| parse_std_align(segs, &align(), 1).unwrap_err();
    for segs in [
        [hit, std_seg(StdLoc::Missing, StdLoc::Missing)],
        [std_seg(StdLoc::Missing, StdLoc::Missing), hit],
    ] {
        assert_eq!(
            refusal(&segs),
            "a Std-seg without locations (NCBI reads past them) is not supported by LOSAT's rpsbproc"
        );
    }
    for only in [int(0, 2, Some(1)), StdLoc::Empty, StdLoc::Whole] {
        assert_eq!(
            refusal(&[hit, std_seg(only, StdLoc::Missing)]),
            "a Std-seg whose only location is not null (NCBI reads past it) is not supported by LOSAT's rpsbproc"
        );
    }
}

// biodata_blast.cpp:178-182 (quoted at `parse_std_align`): `alignedLen = (dst.m_iTo + 1) /
// RF_SIZE` of a minus-strand hit is `int` arithmetic, which overflows for a hit that ends at
// 2^31 - 1 (undefined): its value, 2^31, divided as an `unsigned int` (PD-LOSAT-NCBI-DEFECTS
// exception 11; PS16-Q1 oracle `ntminusmax`), gives 715,827,882 codons, as the hit one base
// shorter does (no overflow). A plus-strand hit with that end parses as it is.
#[test]
fn a_minus_strand_hit_that_ends_at_int_max() {
    let hit = |to: u32, strand: i32| {
        [std_seg(
            int(2_147_483_640, to, Some(strand)),
            int(0, 2, None),
        )]
    };
    let dom = parse_std_align(&hit(0x7fff_ffff, 2), &align(), 1).unwrap();
    assert_eq!(dom.aligned_len(), 715_827_882);
    let dom = parse_std_align(&hit(0x7fff_fffe, 2), &align(), 1).unwrap();
    assert_eq!(dom.aligned_len(), 715_827_882);
    let dom = parse_std_align(&hit(0x7fff_ffff, 1), &align(), 1).unwrap();
    assert_eq!((dom.from, dom.to), (2_147_483_640, i64::from(i32::MAX)));
}

// basedata.cpp:1137-1177 and segset_impl.hpp:281-290 (quoted at `pr2na_ends`): a range of the
// set past `int`'s range (NCBI's set arithmetic overflowed before) and, on the minus strand, a
// residue at `INT_MIN` (`MinusPr2PlusNA`'s `pr - 1`) are computed as their values, which the
// conversions to `unsigned int` take modulo 2^32 (PD-LOSAT-NCBI-DEFECTS exception 11).
#[test]
fn ends_that_overflow_int() {
    let plus = parse_std_align(
        &[std_seg(int(388, 2583, Some(1)), int(0, 731, None))],
        &align(),
        7,
    )
    .unwrap();
    let minus = parse_std_align(
        &[std_seg(int(301, 599, Some(2)), int(0, 99, None))],
        &align(),
        7,
    )
    .unwrap();
    let set = |from: i64, to: i64| SegSet::from_segs(vec![Seg::new(from, to)]);
    let past = i64::from(i32::MAX) + 1;
    // 2^31 converts as `INT_MIN` does: (2^31 * 3 + 1 + 2) modulo 2^32, as an `int`.
    assert_eq!(pr2na_ends(&set(5, past), &plus), Some((16, -2_147_483_645)));
    assert_eq!(
        pr2na_ends(&set(5, past - 1), &plus),
        Some((16, -2_147_483_648))
    );
    assert_eq!(
        pr2na_ends(&set(5, past), &minus),
        Some((-2_147_483_050, 585))
    );
    assert_eq!(
        pr2na_ends(&set(5, past - 1), &minus),
        Some((-2_147_483_047, 585))
    );
    let min = i64::from(i32::MIN);
    // `pr - 1` of `INT_MIN`: -2^31 - 1, 2^31 - 1 modulo 2^32; the residue after it gives 3 less.
    assert_eq!(
        pr2na_ends(&set(min, 10), &minus),
        Some((568, -2_147_483_048))
    );
    assert_eq!(
        pr2na_ends(&set(min + 1, 10), &minus),
        Some((568, -2_147_483_051))
    );
    // The plus strand has no `pr - 1`: `INT_MIN * 3 + rfidx` wraps.
    assert_eq!(pr2na_ends(&set(min, 10), &plus), Some((-2_147_483_647, 33)));
}

// basedata.cpp:1180-1250 (quoted at `segset::SegSet::map_through`): `m_vecSStarts[idx] +
// m_vecLens[idx] - 1` is an `unsigned int` sum stored in the `int` `end`. A segment from model
// position 2^31 - 3 of 10 residues ends at -2^31 + 6 for NCBI, below the site's range (2^31 - 8
// to 2^31 - 1), which is then not mapped; one from 2^31 - 14 ends at 2^31 - 5 and maps the part
// of the range up to it.
#[test]
fn mapped_segment_ends_wrap_as_unsigned_int() {
    let map = |s_start: i32| {
        let dom = parse_align(&dense_seg(&[0, s_start], &[10]), &align(), 1);
        let mut set = SegSet::new();
        set.add(2_147_483_640, 2_147_483_647);
        map_seg_set(&mut set, &dom);
        set.segs()
            .iter()
            .map(|seg| (seg.from, seg.to))
            .collect::<Vec<_>>()
    };
    assert_eq!(map(2_147_483_645), vec![]);
    assert_eq!(map(2_147_483_634), vec![(6, 9)]);
}

// basedata.cpp:1252-1300 (quoted at `translated_pos_map`): a minus-strand hit's positions move
// by `aaLen - alignedLen`; an aligned length beyond the query is NCBI's exception.
#[test]
fn positions_on_a_minus_strand_translation() {
    let segs = [std_seg(int(301, 599, Some(2)), int(0, 99, None))];
    let dom = parse_std_align(&segs, &align(), 7).unwrap();
    let mut set = SegSet::new();
    set.add(10, 12);
    map_seg_set(&mut set, &dom);
    // aaLen = (5000 - 1) / 3 = 1666, alignedLen 200: offset 1466.
    let positions = translated_pos_map(&set, &dom, 5000).unwrap();
    assert_eq!(
        positions
            .iter()
            .map(|p| (p.curr, p.ori))
            .collect::<Vec<_>>(),
        vec![(1476, 10), (1477, 11), (1478, 12)]
    );
    assert_eq!(
        translated_pos_map(&set, &dom, 500),
        Err("Invalid protein length 166: shorter than aligned range 200".to_string())
    );
    // Pr2NaConvert: na length 200 * 3 + 1 = 601; from = 601 - 11 * 3 - 3 - 1 - 3 + 1,
    // to = 601 - 9 * 3 - 3 - 1.
    assert_eq!(pr2na_ends(&set, &dom), Some((562, 570)));
}

// Pr2NaConvert of a plus-strand hit: `pr * 3 + rfidx`, the end `+ 2`.
#[test]
fn ends_on_the_plus_strand() {
    let segs = [std_seg(int(388, 2583, Some(1)), int(0, 731, None))];
    let dom = parse_std_align(&segs, &align(), 7).unwrap();
    let mut set = SegSet::new();
    set.add(5, 9);
    map_seg_set(&mut set, &dom);
    // m_starts 129, rfidx 1: 134 * 3 + 1 = 403, 138 * 3 + 1 + 2 = 417.
    assert_eq!(pr2na_ends(&set, &dom), Some((403, 417)));
    let positions = translated_pos_map(&set, &dom, 5000).unwrap();
    assert_eq!(positions[0].curr, 134);
}
