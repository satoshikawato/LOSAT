use super::*;

// Probes E01-E08, F01-F04 of line DI-8 (`search_strategy.cpp:458-469`).
#[test]
fn the_archive_parameter_follows_the_list_option() {
    use crate::report::archive::options::ParamValue;
    let filter = |taxids: bool, negative: bool| DbFilterArgs {
        taxids: (taxids && !negative).then(|| "1".to_string()),
        negative_taxids: (taxids && negative).then(|| "1".to_string()),
        seqidlist: (!taxids && !negative).then(|| "f".to_string()),
        negative_seqidlist: (!taxids && negative).then(|| "f".to_string()),
        ..DbFilterArgs::default()
    };
    let param = filter(true, false)
        .archive_parameter(&[9606, 10090])
        .unwrap();
    assert_eq!(param.name, "TaxidList");
    assert_eq!(param.value, ParamValue::IntegerList(vec![9606, 10090]));
    // A missing taxonomy file: the list object of NCBI has no id.
    let param = filter(true, false).archive_parameter(&[]).unwrap();
    assert_eq!(param.value, ParamValue::IntegerList(Vec::new()));
    let param = filter(false, true).archive_parameter(&[]).unwrap();
    assert_eq!(param.name, "NegativeTaxidList");
    assert_eq!(param.value, ParamValue::IntegerList(Vec::new()));
    assert!(DbFilterArgs::default().archive_parameter(&[]).is_none());
}

#[test]
fn lines_end_as_ncbis_stream_reader_ends_them() {
    let lines = |text: &str| -> Vec<String> {
        stream_lines(text.as_bytes())
            .into_iter()
            .map(|line| String::from_utf8(line).unwrap())
            .collect()
    };
    assert_eq!(lines("9606\n10090\n"), ["9606", "10090"]);
    assert_eq!(
        lines("9606\r\n\r\n  10090  \r\n\n"),
        ["9606", "", "  10090  ", ""]
    );
    assert_eq!(lines("9606"), ["9606"]);
    assert_eq!(lines("1\r2"), ["1", "2"]);
    assert!(lines("").is_empty());
}
