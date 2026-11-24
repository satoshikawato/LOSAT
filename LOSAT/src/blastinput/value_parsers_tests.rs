// NCBI 2.17.0: `-threshold`'s constraint reads the value with the default flags (`blastp
// -threshold 1e` and `12345678901234567890e`: "Illegal value, expected >=0", exit 1; `0x10`
// searches), `-perc_identity`'s with `fDecimalPosixOrLocal` (`1e` and `0x10` search, `0x1p7`
// is illegal): NCBI's argument reader checks them (`cli::reader`), the conversion reads the number.
#[test]
fn each_constraint_reads_the_value_with_its_own_conversion() {
    let read = |program: &str, name: &str, value: &str| {
        let argv = [
            "losat", program, "-query", "q", "-subject", "s", name, value,
        ];
        crate::cli::try_parse_from(argv)
    };
    let illegal = |program: &str, name: &str, value: &str| match read(program, name, value) {
        Err(crate::cli::ArgvError::Report(stderr)) => {
            String::from_utf8_lossy(&stderr).contains("Illegal value")
        }
        _ => false,
    };
    // The values that the conversion reads after the reader passed them.
    let threshold = |value: &str| match read("blastp", "-threshold", value) {
        Ok(crate::cli::Cli {
            command: crate::cli::Commands::Blastp(args),
        }) => args.threshold.unwrap(),
        other => panic!("{value}: {other:?}"),
    };
    assert_eq!(threshold("0x10"), 16.0);
    assert_eq!(threshold("11"), 11.0);
    assert!(threshold("1e400").is_infinite());
    for value in [
        "1e",
        "2.5E",
        "12345678901234567890e",
        "-nan",
        "+nan(5)",
        "-1",
        "-0x1",
    ] {
        assert!(illegal("blastp", "-threshold", value), "{value}");
    }
    assert!(!illegal("blastp", "-threshold", "1x"));
    assert!(read("blastp", "-threshold", "1x").is_err());
    let percentage = |value: &str| match read("blastn", "-perc_identity", value) {
        Ok(crate::cli::Cli {
            command: crate::cli::Commands::Blastn(args),
        }) => args.percent_identity.unwrap(),
        other => panic!("{value}: {other:?}"),
    };
    assert_eq!(percentage("1e"), 1.0);
    assert_eq!(percentage("0x10"), 16.0);
    match read("blastp", "-qcov_hsp_perc", "99.5E") {
        Ok(crate::cli::Cli {
            command: crate::cli::Commands::Blastp(args),
        }) => assert_eq!(args.qcov_hsp_perc, Some(99.5)),
        other => panic!("-qcov_hsp_perc 99.5E: {other:?}"),
    }
    for value in ["0x1p7", "1e400", "-nan", "-0.5e"] {
        assert!(illegal("blastn", "-perc_identity", value), "{value}");
    }
}
