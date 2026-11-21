use super::*;
use crate::algorithm::blastx::input::parse_fasta;

// No NCBI counterpart: LOSAT's encoding of `-subject` records (`encoded_with_sentinel` is the
// reference).
/// EXPERIMENT (LOSAT_X_BXUNGAPSENT): a FASTA subject read with both sentinels is the record's
/// encoding with the byte 0 added before it, for records of 0, 1 and more letters, lower case
/// and letters outside the alphabet included.
#[test]
fn fasta_subjects_with_both_sentinels_add_the_byte_before() {
    let records =
        parse_fasta(b">a\nMKV*LLAx\n>b\nW\n>c\n\n>d\nbjouzBJOUZ-\n", true, false).unwrap();
    let subjects = Subjects::Fasta(&records);
    for oid in 0..records.len() {
        let both = subjects.x_encoded_with_sentinels(oid);
        assert_eq!(
            *both,
            [&[0u8][..], &subjects.encoded_with_sentinel(oid)].concat(),
            "{oid}"
        );
    }
}

// NCBI reference (598d8ae6): `search_residue` (c++/src/algo/blast/api/blast_setup_cxx.cpp:894-899).
/// A `-subject` record's pyrrolysine (O or o) is searched as X, through every encoding; the
/// record, whose letters the reports show, keeps the O.
#[test]
fn a_fasta_subject_pyrrolysine_is_searched_as_x() {
    let records = parse_fasta(b">a\nSOSoXU\n", true, false).unwrap();
    let subjects = Subjects::Fasta(&records);
    let [s, x, u] = [b'S', b'X', b'U'].map(aa_char_to_ncbistdaa);
    assert_eq!(*subjects.encoded(0), [s, x, s, x, x, u]);
    assert_eq!(*subjects.encoded_with_sentinel(0), [s, x, s, x, x, u, 0]);
    assert_eq!(
        *subjects.x_encoded_with_sentinels(0),
        [0, s, x, s, x, x, u, 0]
    );
    assert_eq!(records[0].sequence, b"SOSOXU");
}
