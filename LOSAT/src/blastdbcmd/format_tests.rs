use super::*;

#[test]
fn the_spec_is_split_as_ncbi_s_loop_does() {
    let spec = FormatSpec::parse("%a%%%l|%%", Vec::new()).unwrap();
    assert_eq!(spec.types, b"al");
    assert_eq!(
        spec.separators,
        [b"".to_vec(), b"%".to_vec(), b"|%".to_vec()]
    );
    assert!(spec.get_defline);
    let spec = FormatSpec::parse("%o %l %h", Vec::new()).unwrap();
    assert!(!spec.get_defline && spec.need_hash && !spec.need_seq);
    // A `%` at the end reads the terminating NUL.
    assert_eq!(FormatSpec::parse("a%", Vec::new()).unwrap().types, [0]);
    assert_eq!(
        FormatSpec::parse("%", Vec::new()).unwrap().separators.len(),
        2
    );
    assert_eq!(
        FormatSpec::parse("no spec", Vec::new()).unwrap_err(),
        "Invalid format specification"
    );
    assert_eq!(
        FormatSpec::parse("%%", Vec::new()).unwrap_err(),
        "Invalid format specification"
    );
    assert_eq!(
        FormatSpec::parse("", Vec::new()).unwrap_err(),
        "Invalid format specification"
    );
}
