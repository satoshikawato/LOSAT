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
