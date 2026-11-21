use std::path::Path;

use super::*;
use crate::blastinput::db_mask::MaskKind;
use crate::seqdb::{MaskRanges, SeqDb, SeqType};

fn mask_prot_v4() -> SeqDb {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db/mask_prot_v4");
    SeqDb::open(path.to_str().unwrap(), SeqType::Protein).unwrap()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1518-1539
// The mask letter covers each mask inside the sequence; a mask past the end is cut there, and
// the masks from the first that starts at or past the end on are not read.
#[test]
fn hard_masked_residues_are_x() {
    let mut residues = vec![1u8; 10];
    mask_protein_letters(&mut residues, &[(0, 0), (2, 4), (6, 20)]);
    assert_eq!(residues, [1, 1, 21, 21, 1, 1, 21, 21, 21, 21]);
    let mut residues = vec![1u8; 4];
    mask_protein_letters(&mut residues, &[(1, 2), (4, 6), (0, 1)]);
    assert_eq!(residues, [1, 21, 1, 1]);
    let mut residues = vec![1u8; 4];
    mask_protein_letters(&mut residues, &[]);
    assert_eq!(residues, [1; 4]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1586-1601
// The closed range `[offset, end]` of an HSP meets a mask that starts at its end (one past its
// last letter); an empty mask holds no letter.
#[test]
fn shown_masks_meet_the_closed_range_of_an_hsp() {
    let masks = [(10, 20), (21, 30), (5, 5), (0, 2), (40, 50)];
    assert_eq!(
        shown_subject_masks(&masks, &[(2, 10)]),
        [(10, 20)],
        "a mask starting at the HSP's end"
    );
    assert_eq!(
        shown_subject_masks(&masks, &[(1, 10)]),
        [(10, 20), (0, 2)],
        "in the order stored"
    );
    assert_eq!(
        shown_subject_masks(&masks, &[(30, 35)]),
        [(21, 30)].repeat(0)
    );
    assert_eq!(
        shown_subject_masks(&masks, &[(29, 35), (45, 46)]),
        [(21, 30), (40, 50)]
    );
    assert!(shown_subject_masks(&masks, &[]).is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-311
// A protein's hard chunks start at their range's left; the soft ranges are the complement in
// the order stored.
#[test]
fn plans_of_soft_and_hard_masks() {
    let mut plan = SubjectPlan::default();
    assert_eq!(
        plan.soft(&[(10, 20), (5, 8)], 30),
        [(0, 10), (20, 5), (8, 30)]
    );
    plan.start_hard(&[(0, 4), (10, 20)], 30);
    let mut chunks = Vec::new();
    while let Some(chunk) = plan.next_chunk() {
        chunks.push((chunk.offset, chunk.length, chunk.residual, chunk.full_range));
    }
    assert_eq!(
        chunks,
        [(0, 0, 0, false), (4, 6, 0, false), (20, 10, 0, false)]
    );
    plan.start_hard(&[(10, 30)], 30);
    let chunk = plan.next_chunk().unwrap();
    assert!(chunk.full_range && chunk.length == 10, "{chunk:?}");
    assert!(plan.next_chunk().is_none());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-500 and 762-766 (cited at
// `SubjectPlan::next_searched_chunk`): a chunk shorter than the word gives nothing in either
// scan, so the searched chunks are those that hold a word.
#[test]
fn chunks_shorter_than_the_word_are_passed_over() {
    let searched = |masks: &[(u32, u32)], length: usize, word: usize| {
        let mut plan = SubjectPlan::default();
        plan.start_hard(masks, length);
        let mut chunks = Vec::new();
        while let Some(chunk) = plan.next_searched_chunk(word) {
            chunks.push((chunk.offset, chunk.length));
        }
        chunks
    };
    // The hand masks of p5 (`dbmaskp/edge_tiny.asnt`): masks of one letter at 101 and 104 leave
    // the two letters 102 and 103, passed over at a word of 3 and searched at a word of 2.
    let tiny = [(101, 102), (104, 105), (200, 202)];
    assert_eq!(searched(&tiny, 250, 3), [(0, 101), (105, 95), (202, 48)]);
    assert_eq!(
        searched(&tiny, 250, 2),
        [(0, 101), (102, 2), (105, 95), (202, 48)]
    );
    // Touching masks leave a chunk of no letters, a mask at the start one at 0, and a chunk of
    // exactly the word's letters is searched (`dbmaskp/edge_short.asnt`).
    assert_eq!(searched(&[(30, 60), (60, 90)], 150, 3), [(0, 30), (90, 60)]);
    assert_eq!(searched(&[(0, 4), (10, 20)], 30, 3), [(4, 6), (20, 10)]);
    assert_eq!(searched(&[(1, 10)], 150, 3), [(10, 140)]);
    assert_eq!(
        searched(&[(30, 60), (63, 90)], 150, 3),
        [(0, 30), (60, 3), (90, 60)]
    );
    assert_eq!(searched(&[(30, 60), (62, 90)], 150, 3), [(0, 30), (90, 60)]);
    // Every chunk short: nothing is searched.
    assert!(searched(&[(2, 5), (7, 9)], 10, 3).is_empty());
}

// The reader gives each sequence's stored masks as `SeqDb::mask_ranges_at` does, across the
// windows it reads, and its search's kind and algorithm come from the resolved option.
#[test]
fn the_reader_gives_the_stored_masks_of_each_sequence() {
    let db = mask_prot_v4();
    let ids = db.mask_algorithms().unwrap().ids();
    assert!(ids.contains(&21), "{ids:?}");
    for window in [1, 2, 4096] {
        let mut reader = SubjectMasks::with_window(window);
        let mut with_masks = 0;
        for ordinal in (0..db.num_seqs()).chain([1, 0]) {
            let expected = match db.mask_ranges_at(db.locate_oid(ordinal).unwrap(), 21) {
                Ok(MaskRanges::Ranges(ranges)) => ranges,
                other => panic!("{ordinal}: {other:?}"),
            };
            let found = reader
                .of(
                    &db,
                    ordinal,
                    21,
                    crate::algorithm::common::database::search_read_error_of,
                    "BLASTP",
                )
                .unwrap()
                .to_vec();
            assert_eq!(found, expected, "window {window}, sequence {ordinal}");
            with_masks += usize::from(!found.is_empty());
        }
        assert!(with_masks > 0);
    }
    let resolved = ResolvedMask {
        kind: MaskKind::Hard,
        id: 21,
        key: None,
        shown: Vec::new(),
    };
    let mask = SearchMask::of(Some(&resolved)).unwrap();
    assert_eq!(mask.algo_id, 21);
    assert!(mask.hard());
    assert_eq!(SearchMask::of(None), None);
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1518-1539, 1672
// A nucleotide's copied fetch writes ncbi4na B over every masked base after its ambiguities, in
// the order of the protein letters' masks (masks of no letter at 0 skipped, one past the end cut,
// a later mask out of order applied, the masks from the first that starts at or past the end on
// not read); the packed bases stay those the preliminary stage reads.
#[test]
fn hard_masked_bases_are_b_over_the_ambiguities() {
    let store =
        crate::sequence::packed_subject::PackedNucleotide::from_subject_letters(b"ACGTNACGTRACGT")
            .unwrap();
    let letters = |store: &crate::sequence::packed_subject::PackedNucleotide| {
        let mut out = Vec::new();
        store.ascii_into(&mut out);
        String::from_utf8(out).unwrap()
    };
    assert_eq!(letters(&store), "ACGTNACGTRACGT");
    let masked = mask_nucleotide_store(&store, &[(0, 0), (1, 3), (4, 6), (12, 20), (0, 1)]);
    assert_eq!(letters(&masked), "BBBTBBCGTRACBB");
    assert_eq!(masked.packed(), store.packed());
    assert_eq!(masked.len(), store.len());
    // A mask that starts at the end stops the others; no masks change nothing.
    let masked = mask_nucleotide_store(&store, &[(14, 15), (0, 2)]);
    assert_eq!(letters(&masked), "ACGTNACGTRACGT");
    assert_eq!(
        letters(&mask_nucleotide_store(&store, &[])),
        "ACGTNACGTRACGT"
    );
    assert_eq!(HARD_MASK_NUCLEOTIDE_CODE, 14);
}

/// A pairwise row of query `q_idx` with subject rank `s_idx` over the subject's letters
/// `s_start..=s_end` (1-based).
fn row(
    q_idx: u32,
    s_idx: u32,
    s_start: usize,
    s_end: usize,
) -> crate::report::pairwise::PairwiseHit {
    crate::report::pairwise::PairwiseHit::from(crate::common::Hit {
        identity: 100.0,
        length: s_end + 1 - s_start,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: s_end + 1 - s_start,
        s_start,
        s_end,
        e_value: 1e-10,
        bit_score: 50.0,
        num_ident: s_end + 1 - s_start,
        query_frame: 0,
        query_length: 400,
        q_idx,
        s_idx,
        raw_score: 100,
        sort_query_offset: 0,
        sort_query_end: s_end + 1 - s_start,
        sort_subject_offset: s_start - 1,
        sort_subject_end: s_end,
        has_sort_offsets: true,
        gap_info: None,
        num_positives: s_end + 1 - s_start,
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1046-1052 and 1247-1250
// (cited at `shown_db_subject_masks`)
// Two sequences of one Seq-id (ranks 0 and 1, from two databases with parsed Seq-ids): both rows
// of a query show the masks of the later in the rows' order that meet its own HSPs; a sequence
// whose masks meet none of its HSPs replaces nothing; another Seq-id keeps its own; each query
// keeps its own map.
#[test]
fn rows_show_the_masks_of_their_seq_id() {
    let stored: [Vec<(u32, u32)>; 3] = [vec![(150, 190)], vec![(0, 30), (300, 320)], vec![(5, 9)]];
    let keys = |rank: u32| match rank {
        0 | 1 => (b"p5".to_vec(), b"p5".to_vec()),
        _ => (b"p6".to_vec(), b"p6".to_vec()),
    };
    let hits = [
        row(0, 1, 21, 380),
        row(0, 0, 21, 380),
        row(0, 2, 1, 50),
        row(1, 0, 140, 200),
        row(1, 1, 100, 120),
    ];
    let shown = shown_db_subject_masks(
        &hits,
        |rank| stored.get(rank as usize).map(Vec::as_slice),
        keys,
    );
    // Query 0: rank 1's [0, 30) meets its HSP (closed [20, 380]), [300, 320) too; rank 0's
    // [150, 190) comes later and replaces them for both rows.
    assert_eq!(shown.get(&(0, 0)), Some(&vec![(150, 189)]));
    assert_eq!(shown.get(&(0, 1)), Some(&vec![(150, 189)]));
    assert_eq!(shown.get(&(0, 2)), Some(&vec![(5, 8)]));
    // Query 1: rank 0's mask meets its HSP; rank 1's masks meet none of [99, 120], so rank 0's
    // stay for both rows.
    assert_eq!(shown.get(&(1, 0)), Some(&vec![(150, 189)]));
    assert_eq!(shown.get(&(1, 1)), Some(&vec![(150, 189)]));
    // Ranks without a shared Seq-id: each its own (the BL_ORD_ID key of each OID).
    let own = |rank: u32| (vec![rank as u8], vec![rank as u8]);
    let shown = shown_db_subject_masks(
        &hits[..2],
        |rank| stored.get(rank as usize).map(Vec::as_slice),
        own,
    );
    assert_eq!(shown.get(&(0, 1)), Some(&vec![(0, 29), (300, 319)]));
    assert_eq!(shown.get(&(0, 0)), Some(&vec![(150, 189)]));
}
