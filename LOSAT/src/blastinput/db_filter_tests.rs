use super::*;

// Probes E01-E08, F01-F04 of line DI-8 (`search_strategy.cpp:444-498`) and D6 GI/IPG lists.
#[test]
fn the_archive_parameter_follows_the_list_option() {
    use crate::report::archive::options::{Param, ParamValue};
    let filter = |taxids: bool, negative: bool| DbFilterArgs {
        taxids: (taxids && !negative).then(|| "1".to_string()),
        negative_taxids: (taxids && negative).then(|| "1".to_string()),
        seqidlist: (!taxids && !negative).then(|| "f".to_string()),
        negative_seqidlist: (!taxids && negative).then(|| "f".to_string()),
        ..DbFilterArgs::default()
    };
    assert_eq!(
        filter(true, false).archive_parameters(&[], &[9606, 10090]),
        vec![Param {
            name: "TaxidList",
            value: ParamValue::IntegerList(vec![9606, 10090]),
        }]
    );
    // A missing taxonomy file: the list object of NCBI has no id.
    assert_eq!(
        filter(true, false).archive_parameters(&[], &[]),
        vec![Param {
            name: "TaxidList",
            value: ParamValue::IntegerList(Vec::new()),
        }]
    );
    assert_eq!(
        filter(false, true).archive_parameters(&[], &[]),
        vec![Param {
            name: "NegativeTaxidList",
            value: ParamValue::IntegerList(Vec::new()),
        }]
    );
    assert!(DbFilterArgs::default()
        .archive_parameters(&[], &[])
        .is_empty());

    // Positive and negative GI lists write GiList/NegativeGiList + empty TaxidList/NegativeTaxidList,
    // and a 0-GI file writes only the empty TaxidList/NegativeTaxidList.
    let pos_gi = DbFilterArgs {
        gilist: Some("g.gil".to_string()),
        ..DbFilterArgs::default()
    };
    assert_eq!(
        pos_gi.archive_parameters(&[101, 102], &[]),
        vec![
            Param {
                name: "GiList",
                value: ParamValue::BigIntegerList(vec![101, 102]),
            },
            Param {
                name: "TaxidList",
                value: ParamValue::IntegerList(Vec::new()),
            },
        ]
    );
    assert_eq!(
        pos_gi.archive_parameters(&[], &[]),
        vec![Param {
            name: "TaxidList",
            value: ParamValue::IntegerList(Vec::new()),
        }]
    );
    let neg_gi = DbFilterArgs {
        negative_gilist: Some("g.gil".to_string()),
        ..DbFilterArgs::default()
    };
    assert_eq!(
        neg_gi.archive_parameters(&[102, 101], &[]),
        vec![
            Param {
                name: "NegativeGiList",
                value: ParamValue::BigIntegerList(vec![102, 101]),
            },
            Param {
                name: "NegativeTaxidList",
                value: ParamValue::IntegerList(Vec::new()),
            },
        ]
    );

    // IPG lists write empty TaxidList / NegativeTaxidList; other lists take precedence over IPG.
    let pos_ipg = DbFilterArgs {
        ipglist: Some("p.ipg".to_string()),
        negative_ipglist: Some("n.ipg".to_string()),
        ..DbFilterArgs::default()
    };
    assert_eq!(
        pos_ipg.archive_parameters(&[], &[]),
        vec![Param {
            name: "TaxidList",
            value: ParamValue::IntegerList(Vec::new()),
        }]
    );
    let neg_ipg = DbFilterArgs {
        negative_ipglist: Some("n.ipg".to_string()),
        ..DbFilterArgs::default()
    };
    assert_eq!(
        neg_ipg.archive_parameters(&[], &[]),
        vec![Param {
            name: "NegativeTaxidList",
            value: ParamValue::IntegerList(Vec::new()),
        }]
    );
    let gi_and_ipg = DbFilterArgs {
        gilist: Some("g.gil".to_string()),
        ipglist: Some("p.ipg".to_string()),
        ..DbFilterArgs::default()
    };
    assert_eq!(
        gi_and_ipg.archive_parameters(&[101], &[]),
        vec![
            Param {
                name: "GiList",
                value: ParamValue::BigIntegerList(vec![101]),
            },
            Param {
                name: "TaxidList",
                value: ParamValue::IntegerList(Vec::new()),
            },
        ]
    );
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
