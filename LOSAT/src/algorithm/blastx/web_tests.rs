use super::*;

// ABI v1 is frozen (plan TD-1): its BLASTX `-outfmt` keeps the values and the parser
// errors of the CLI parser before BLASTX read `-outfmt` as NCBI (dev 0797cc40, whose
// `losat blastx -outfmt 5` wrote these lines with exit status 2).
#[test]
fn web_v1_keeps_its_outfmt_values_and_errors() {
    let error = |outfmt: &str| {
        let error = parse_args(&["-outfmt", outfmt]).unwrap_err();
        crate::cli::render_message(error.downcast_ref::<clap::Error>().unwrap())
    };
    for (outfmt, reason) in [
        ("5", "unsupported BLASTX outfmt: use 0, 6 or 7"),
        ("abc", "unsupported BLASTX outfmt: use 0, 6 or 7"),
        ("13", "unsupported BLASTX outfmt: use 0, 6 or 7"),
        ("0 std", "custom fields require outfmt 6/7"),
        ("6 delim=", "unsupported BLASTX outfmt field 'delim='"),
        ("7 qaccver foo", "unsupported BLASTX outfmt field 'foo'"),
    ] {
        assert_eq!(
            error(outfmt),
            format!(
                "error: invalid value '{outfmt}' for '-outfmt <OUTFMT>': {reason}\n\nFor more information, try '--help'.\n"
            )
        );
    }
    // White space between the words, read as ABI v1 read it.
    for (outfmt, number, fields) in [
        ("6\tqseqid", 6, &["qseqid"][..]),
        (" 7 ", 7, &[]),
        ("0", 0, &[]),
    ] {
        let args = parse_args(&["-outfmt", outfmt]).unwrap();
        let options = v1_resolve(&args).unwrap();
        assert_eq!(options.outfmt, number, "{outfmt:?}");
        if !fields.is_empty() {
            let spec = options.tabular.expect("a field list");
            let names: Vec<&str> = spec.fields.iter().map(|f| f.name()).collect();
            assert_eq!(names, fields, "{outfmt:?}");
        }
    }
}
