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

/// The counts of an ascii counts file whose unit size line is `unit_size` (read as a `Uint1`).
fn ascii_counts(unit_size: &str) -> UnitCounts {
    let text =
        format!("{unit_size}\n0 5\n3 7\n\n>t_threshold 6\n>t_extend 5\n>t_low 1\n>t_high 100\n");
    let params = super::super::MaskerParams::default();
    UnitCounts::read(&mut io::Cursor::new(text.into_bytes()), &params).unwrap()
}

/// `convert_counts` of `istat` into `format`, then `finalize`.
fn convert(istat: &UnitCounts, format: CountsFormat) -> Result<Vec<u8>, FinalizeError> {
    let mut writer = CountsWriter::new(format, b"");
    convert_counts(istat, Vec::new(), &mut writer, &mut |_| {}).unwrap();
    let mut out = Vec::new();
    writer.finalize(&mut out, &mut |_| {}).map(|()| out)
}

#[test]
fn optimized_formats_of_units_of_0_letters_are_refused() {
    let error = optimize(0, &[0], &[1], 1, None).unwrap_err();
    assert!(matches!(error, FinalizeError::Refused(_)));
    for size in ["0", "256"] {
        let istat = ascii_counts(size);
        assert_eq!(istat.unit_size(), 0);
        // NCBI writes the ascii and binary formats of them: the unit 0 alone.
        let ascii = convert(&istat, CountsFormat::Ascii).unwrap();
        assert!(ascii.starts_with(b"##windowmasker-statistics-format-version:ascii"));
        assert!(ascii.windows(5).any(|line| line == b"\n0\n0 "));
        assert!(convert(&istat, CountsFormat::Binary).is_ok());
        for format in [
            CountsFormat::OAscii {
                size_requested: 512,
            },
            CountsFormat::OBinary {
                size_requested: 512,
            },
        ] {
            assert!(
                matches!(convert(&istat, format), Err(FinalizeError::Refused(_))),
                "{format:?}"
            );
        }
    }
}

#[test]
fn convert_of_units_of_17_letters_is_refused() {
    let istat = ascii_counts("17");
    let mut writer = CountsWriter::new(CountsFormat::Ascii, b"");
    let error = convert_counts(&istat, Vec::new(), &mut writer, &mut |_| {}).unwrap_err();
    assert!(matches!(error, WinMaskError::UnitSize { unit: 17 }));
}
