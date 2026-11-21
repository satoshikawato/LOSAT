use super::*;
use crate::blastinput::fasta_reader::{read_subjects, FastaInputSource, ReaderConfig};
use crate::blastinput::seq_range::{cut_subjects, SequenceRange};

// No NCBI counterpart: when the letters of a subject are packed for its plan. The test holds
// with the switch on (the default) and off (`LOSAT_X_TBNSTORE=0` in the test's environment).
#[test]
fn a_plan_packs_its_letters_unless_it_reads_their_case_or_cannot() {
    let on = x_tbn_store();
    let packed = |letters: &[u8], mask_lowercase: bool| {
        x_plan_store(letters, mask_lowercase).map(|store| {
            let mut back = Vec::new();
            store.ascii_into(&mut back);
            back
        })
    };
    // every IUPAC letter, `-` and lower case: the store gives the letters back
    let letters = b"ACGTacgtNRYKMSWBDHVnrykmswbdhv-ACGT".to_vec();
    assert_eq!(packed(&letters, false), on.then(|| letters.clone()));
    // `-lcase_masking` reads the case: letters with lower case stay letters
    assert_eq!(packed(&letters, true), None);
    let upper = b"ACGTNRYKMSWBDHV".to_vec();
    assert_eq!(packed(&upper, true), on.then(|| upper.clone()));
    // a byte without an ncbi4na code (only an API caller can give one)
    assert_eq!(packed(b"ACGU", false), None);
    assert_eq!(packed(b"", false), on.then(Vec::new));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-431
// NCBI packs the subject it searches, the interval of `-subject_loc`, and seeds the resolution
// of its ambiguous bases with the interval's length. With a range the reader packs the records
// `cut_subjects` makes (at one thread and on the packing thread): their letters, deflines,
// placements and resolved bases; the whole record packed and read in the interval resolves
// the ambiguous bases otherwise.
#[test]
fn a_range_packs_the_records_cut_to_it() {
    let mut state = 5u64;
    let mut fasta = Vec::new();
    for (index, length) in [3_000usize, 2_500, 4_000].into_iter().enumerate() {
        fasta.extend_from_slice(format!(">s{index} title {index}\n").as_bytes());
        for _ in 0..length {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let draw = (state >> 33) as usize;
            fasta.push(if draw.is_multiple_of(8) {
                b"RYKMSWBDHVNrykmswbdhvn"[draw / 8 % 22]
            } else {
                b"ACGTacgt"[draw / 8 % 8]
            });
        }
        fasta.push(b'\n');
    }
    let range = SequenceRange {
        from: 100,
        to: 2_800,
    };
    let config = || ReaderConfig::subject("TBLASTN", false, false);
    let whole = read_subjects(
        &mut FastaInputSource::from_bytes(&fasta, config()),
        Some(&range),
        &mut |_| Ok(()),
    )
    .unwrap();
    let (cut, expected_placements) = cut_subjects(&whole, Some(&range)).unwrap().unwrap();
    let mut resolved_otherwise = 0;
    for threads in [1, 2] {
        let (records, stores, placements) = x_read_subject_stores(
            &mut FastaInputSource::from_bytes(&fasta, config()),
            Some(&range),
            &mut |_| Ok(()),
            threads,
        )
        .unwrap();
        let stores = stores.expect("every record packed");
        assert_eq!(placements.0, expected_placements.0, "{threads} threads");
        assert_eq!(records.len(), cut.len());
        for (index, ((record, store), cut)) in records.iter().zip(&stores).zip(&cut).enumerate() {
            assert_eq!(
                (&record.local_id, &record.title, record.sequence.is_empty()),
                (&cut.local_id, &cut.title, true)
            );
            let (mut letters, mut resolved, mut expected) = (Vec::new(), Vec::new(), Vec::new());
            store.ascii_into(&mut letters);
            assert_eq!(letters, cut.sequence, "{threads} threads, record {index}");
            store.resolved_letters_into(&mut resolved);
            PackedNucleotide::from_subject_letters(&cut.sequence)
                .unwrap()
                .resolved_letters_into(&mut expected);
            assert_eq!(resolved, expected, "{threads} threads, record {index}");
            let mut whole_resolved = Vec::new();
            PackedNucleotide::from_subject_letters(&whole[index].sequence)
                .unwrap()
                .resolved_letters_into(&mut whole_resolved);
            let offset = placements.offset(index);
            if whole_resolved[offset..offset + resolved.len()] != resolved[..] {
                resolved_otherwise += 1;
            }
        }
    }
    assert!(resolved_otherwise > 0, "the seed is the interval's length");
}
