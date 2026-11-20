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
