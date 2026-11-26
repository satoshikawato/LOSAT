use super::*;

// ==========================================================================
// E-value formatting tests (NCBI parity)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:965-988
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1328-1337
// ==========================================================================
// ==========================================================================

#[test]
fn test_format_evalue_ncbi_zero() {
    // Values at or below 0.0 should return "0.0"
    assert_eq!(format_evalue_ncbi_tabular(0.0), "0.0");
}

#[test]
fn test_format_evalue_ncbi_very_small() {
    // < 1e-180: "0.0"
    assert_eq!(format_evalue_ncbi_tabular(1e-200), "0.0");
    assert_eq!(format_evalue_ncbi_tabular(1e-181), "0.0");
}

#[test]
fn test_format_evalue_ncbi_small_scientific() {
    // [1e-180, 0.0009): tabular override uses "%.2e"
    assert_eq!(format_evalue_ncbi_tabular(1e-100), "1.00e-100");
    assert_eq!(format_evalue_ncbi_tabular(1e-150), "1.00e-150");
}

#[test]
fn test_format_evalue_ncbi_medium_scientific() {
    // [1e-180, 0.0009): "%.2e" format (tabular.cpp override)
    assert_eq!(format_evalue_ncbi_tabular(1e-50), "1.00e-50");
    assert_eq!(format_evalue_ncbi_tabular(1.23e-20), "1.23e-20");
    assert_eq!(format_evalue_ncbi_tabular(5e-5), "5.00e-05");
}

#[test]
fn test_format_evalue_ncbi_decimal() {
    // [0.0009, 0.1): "%4.3lf" - 3 decimal places
    assert_eq!(format_evalue_ncbi_tabular(0.001), "0.001");
    assert_eq!(format_evalue_ncbi_tabular(0.005), "0.005");
    assert_eq!(format_evalue_ncbi_tabular(0.099), "0.099");

    // [0.1, 1.0): "%3.2lf" - 2 decimal places
    assert_eq!(format_evalue_ncbi_tabular(0.1), "0.10");
    assert_eq!(format_evalue_ncbi_tabular(0.5), "0.50");
    assert_eq!(format_evalue_ncbi_tabular(0.99), "0.99");

    // [1.0, 10.0): "%2.1lf" - 1 decimal place
    assert_eq!(format_evalue_ncbi_tabular(1.0), "1.0");
    assert_eq!(format_evalue_ncbi_tabular(5.5), "5.5");
    assert_eq!(format_evalue_ncbi_tabular(9.9), "9.9");
}

#[test]
fn test_format_evalue_ncbi_integer() {
    // >= 10.0: "%2.0lf" - no decimal
    assert_eq!(format_evalue_ncbi_tabular(10.0), "10");
    assert_eq!(format_evalue_ncbi_tabular(100.0), "100");
    assert_eq!(format_evalue_ncbi_tabular(1000.0), "1000");
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:965-975
#[test]
fn test_format_evalue_ncbi_base_scientific() {
    // Base GetScoreString formatting uses no decimals for < 0.0009
    assert_eq!(format_evalue_ncbi(1e-100), "1e-100");
    assert_eq!(format_evalue_ncbi(1e-50), "1e-50");
}

// ==========================================================================
// Bit score formatting tests (NCBI parity)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:986-993
// ==========================================================================

#[test]
fn test_format_bitscore_ncbi_small() {
    // <= 99.9: "%4.1lf" - one decimal place
    assert_eq!(format_bitscore_ncbi(0.0), " 0.0");
    assert_eq!(format_bitscore_ncbi(9.8), " 9.8");
    let mut buffer = Vec::new();
    write_bitscore_ncbi(&mut buffer, 9.8).unwrap();
    assert_eq!(buffer, b" 9.8");
    assert_eq!(format_bitscore_ncbi(50.5), "50.5");
    assert_eq!(format_bitscore_ncbi(99.9), "99.9");
}

#[test]
fn test_format_bitscore_ncbi_medium() {
    // > 99.9 and <= 99999: "%3.0ld" - integer (no decimal)
    assert_eq!(format_bitscore_ncbi(100.0), "100");
    // "%3.0ld" pads (long)99.95 = 99 to three characters.
    assert_eq!(format_bitscore_ncbi(99.95), " 99");
    let mut buffer = Vec::new();
    write_bitscore_ncbi(&mut buffer, 99.95).unwrap();
    assert_eq!(buffer, b" 99");
    assert_eq!(format_bitscore_ncbi(185.5), "185");
    assert_eq!(format_bitscore_ncbi(598.9), "598");
    assert_eq!(format_bitscore_ncbi(692.0), "692");
    assert_eq!(format_bitscore_ncbi(99999.0), "99999");
}

#[test]
fn test_format_bitscore_ncbi_large() {
    // > 99999: "%5.3le" - scientific notation
    assert_eq!(format_bitscore_ncbi(100000.0), "1.000e+05");
    assert_eq!(format_bitscore_ncbi(123456.789), "1.235e+05");
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:965-988
#[test]
fn test_exponent_padding_without_generated_fixture() {
    assert_eq!(format_evalue_ncbi_tabular(5e-5), "5.00e-05");
    assert_eq!(format_bitscore_ncbi(100000.0), "1.000e+05");
}

// ==========================================================================
// Legacy tests
// ==========================================================================

#[test]
fn test_format_percent_identity_ncbi_divides_before_scaling() {
    // NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:437-441
    assert_eq!(format_percent_identity_ncbi(179, 320, 3), "55.937");
    assert_eq!(format_percent_identity_ncbi(0, 0, 3), "0.000");
}

// NCBI BLAST+ 2.17.0, `blastn -outfmt 7` (inventory range RP,
// `scratch_RP/out/<case>/blastn_f7.out`): the query's title as the reader keeps it,
// raw, after `# Query: `.
#[test]
fn outfmt7_query_lines_follow_the_ncbi_oracle() {
    use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
    for (defline, line) in [
        (&b"id desc"[..], &b"# Query: id desc\n"[..]),
        (b"id\tdesc", b"# Query: id\n"),
        (b"", b"# Query: \n"),
        (b"id\xe9x desc", b"# Query: id\xe9x desc\n"),
        (
            b"id\xffx desc \xff end",
            b"# Query: id\xffx desc \xff end\n",
        ),
        (
            b"id &amp; &lt;b&gt; \"q\" desc",
            b"# Query: id &amp; &lt;b&gt; \"q\" desc\n",
        ),
        (b"\x01id desc", b"# Query: \x01id desc\n"),
        (b"id desc   ", b"# Query: id desc\n"),
        (b"TPA: id desc", b"# Query: TPA: id desc\n"),
        (b"id desc\x01more", b"# Query: id desc\n"),
    ] {
        let bytes = [&b">"[..], defline, b"\nACGT\n"].concat();
        let mut source =
            FastaInputSource::from_bytes(&bytes, ReaderConfig::query("BLASTN", false, false));
        let title = read_all(&mut source, &mut |_| Ok(())).unwrap()[0]
            .title
            .clone();
        let mut out = Vec::new();
        write_outfmt7_query_line(&mut out, &title).unwrap();
        assert_eq!(out, line, "{defline:?}");
    }
}

// Line PS item PS7c: an infinite bit score (a PSSM whose gapped Lambda is `PLUS-INFINITY`) is
// printed by C's `%5.3le` as `inf` in the field of 5, in both writers.
#[test]
fn an_infinite_bit_score_fills_a_field_of_five() {
    assert_eq!(format_bitscore_ncbi(f64::INFINITY), "  inf");
    let mut out = Vec::new();
    write_bitscore_ncbi(&mut out, f64::INFINITY).unwrap();
    assert_eq!(out, b"  inf");
    assert_eq!(format_bitscore_ncbi(123456.0), "1.235e+05");
}
