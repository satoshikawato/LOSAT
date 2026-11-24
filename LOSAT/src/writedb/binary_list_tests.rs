//! The binary GI lists that NCBI BLAST+ 2.17.0 blastdb_aliastool wrote (probes G03, G13 of row
//! RG-D2a; the fixture rows of `tests/fixtures/blastdb_aliastool` compare whole runs).

use super::gi_list_bytes;
use crate::seqdb::id_list::read_gi_list;

#[test]
fn four_byte_gis_sorted_with_repeats() {
    let bytes = gi_list_bytes(&[105, 101, 105, 5]);
    assert_eq!(
        bytes,
        [
            0xFF, 0xFF, 0xFF, 0xFF, 0, 0, 0, 4, 0, 0, 0, 5, 0, 0, 0, 0x65, 0, 0, 0, 0x69, 0, 0, 0,
            0x69
        ]
    );
    assert_eq!(read_gi_list(&bytes), Ok(vec![5, 101, 105, 105]));
}

#[test]
fn a_negative_or_wide_gi_makes_eight_byte_words() {
    let bytes = gi_list_bytes(&[101, -5, 0, 4_294_967_297]);
    assert_eq!(&bytes[..8], &[0xFF, 0xFF, 0xFF, 0xFE, 0, 0, 0, 4]);
    assert_eq!(&bytes[8..16], &(-5i64).to_be_bytes());
    assert_eq!(&bytes[32..40], &4_294_967_297i64.to_be_bytes());
    assert!(read_gi_list(&bytes).is_err(), "NCBI's GI reader takes -1 lists only");
}

#[test]
fn no_gis() {
    assert_eq!(gi_list_bytes(&[]), [0xFF, 0xFF, 0xFF, 0xFF, 0, 0, 0, 0]);
}
