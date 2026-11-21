use super::*;

#[test]
fn the_names_are_in_order_and_the_values_are_ncbi_s() {
    assert!(ACC_INFO_NAMES.windows(2).all(|pair| pair[0].0 < pair[1].0));
    // Values of the enumerators (NCBI_ACC of Seq_id.hpp).
    assert_eq!(find_acc_info(b"gb_prot"), Some(0x4000_0005));
    assert_eq!(find_acc_info(b"REFSEQ_PROT_PREDICTED"), Some(0x6000_000a));
    assert_eq!(find_acc_info(b"refseq_chromosome"), Some(0x8000_850a));
    assert_eq!(find_acc_info(b"gb_wgsvm_nuc"), Some(0x8900_8305));
    assert_eq!(find_acc_info(b"pdb"), Some(15));
    assert_eq!(find_acc_info(b"unknown"), Some(UNKNOWN));
    assert_eq!(find_acc_info(b"no_such_type"), None);
}
