use super::*;

#[test]
fn no_memory_for_units_of_up_to_three_letters_is_refused() {
    let error = optimize(3, &[0], &[1], 0, None).unwrap_err();
    assert!(matches!(error, FinalizeError::Refused(_)));
    let error = optimize(4, &[0], &[1], 0, None).unwrap_err();
    assert!(matches!(error, FinalizeError::Ncbi(_)));
}

#[test]
fn the_size_of_smem_is_its_low_16_bits() {
    let size = |smem| match CountsFormat::from_sformat("obinary", smem) {
        CountsFormat::OBinary { size_requested } => size_requested,
        _ => unreachable!(),
    };
    assert_eq!(size(65536), 0);
    assert_eq!(size(70000), 4464);
    assert_eq!(size(-1), 65535);
    assert_eq!(size(512), 512);
}
