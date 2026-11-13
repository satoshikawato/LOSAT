use super::*;

fn run_words(line: &[&str], input: &[u8]) -> (i32, String, String) {
    let words: Vec<OsString> = line.iter().map(OsString::from).collect();
    let mut out = Vec::new();
    let mut err = Vec::new();
    let status = run_with(&words, &mut &input[..], &mut out, &mut err);
    (
        status,
        String::from_utf8_lossy(&out).into_owned(),
        String::from_utf8(err).unwrap(),
    )
}

const ALGO: [&str; 4] = ["-masking_algorithm", "dust", "-masking_options", "o"];

fn with_algo<'a>(line: &[&'a str]) -> Vec<&'a str> {
    let mut words = ALGO.to_vec();
    words.extend_from_slice(line);
    words
}

#[test]
fn the_short_usage_is_the_head_of_the_help() {
    let usage = crate::ncbi_app::short_usage(HELP);
    assert_eq!(usage.len(), 414);
    assert!(usage.starts_with("USAGE\n  convert2blastmask [-h] [-help] [-in input_file_name]"));
    assert_eq!(HELP.len(), 1203);
}

#[test]
fn the_version_is_blast_s() {
    let (status, out, err) = run_words(&["-version"], b"");
    assert_eq!((status, err), (0, String::new()));
    assert!(
        out.starts_with(&format!(
            "convert2blastmask: {}\n Package: LOSAT ",
            crate::NCBI_BLAST_VERSION
        )),
        "{out}"
    );
}

// The full versions, as dustmasker's, before the mandatory keys are checked (line RG item E3;
// NCBI BLAST+ 2.17.0 `convert2blastmask -version-full`: `convert2blastmask: 2.17.0+`, rc 0).
#[test]
fn the_full_versions_are_blast_s() {
    use crate::ncbi_version::{print, VersionRequest};
    for (word, request) in [
        ("-version-full", VersionRequest::Full),
        ("-version-full-xml", VersionRequest::FullXml),
        ("-version-full-json", VersionRequest::FullJson),
    ] {
        assert_eq!(
            run_words(&["-out", word], b""),
            (
                0,
                print(request, "convert2blastmask", crate::NCBI_BLAST_VERSION),
                String::new()
            ),
            "{word}"
        );
    }
}

// NCBI 2.17.0: `Argument "masking_algorithm". Mandatory value is missing:  `String'` after
// the usage, exit status 1 (fixture `noalgo`).
#[test]
fn the_mandatory_keys_must_be_given() {
    let (status, out, err) = run_words(&["-masking_options", "o"], b"");
    assert_eq!((status, out), (1, String::new()));
    assert!(err.ends_with(
        "Error: Argument \"masking_algorithm\". Mandatory value is missing:  `String'\n\
         Error:  (CArgException::eNoArg) Argument \"masking_algorithm\". Mandatory value is \
         missing:  `String'\n"
    ));
}

// NCBI 2.17.0 on `>s1 desc / ACGTacgtACGTNNNNnnnnACGT / acgt / >s2 / ACGTACGT / >s3 /
// acgtacgtACGT` (fixture `lc_int`): an interval runs on over the end of a line, and a record
// without a lower case letter has none.
#[test]
fn the_lower_case_letters_are_the_intervals() {
    let input = b">s1 desc\nACGTacgtACGTNNNNnnnnACGT\nacgt\n>s2\nACGTACGT\n>s3\nacgtacgtACGT\n";
    assert_eq!(
        run_words(&with_algo(&["-outfmt", "interval"]), input),
        (
            0,
            "4 - 7\n16 - 19\n24 - 27\n0 - 7\n".to_string(),
            String::new()
        )
    );
}

// The interval line and the `A` in its place: an interval that a lower case letter opened
// closes after the interval line's own (fixture `iv_lc_int`).
#[test]
fn an_interval_line_closes_an_open_interval_after_its_own() {
    assert_eq!(
        run_words(
            &with_algo(&["-outfmt", "interval"]),
            b">e\nacgt\n10 - 20\nacgtACGT\n"
        ),
        (0, "10 - 20\n0 - 3\n5 - 8\n".to_string(), String::new())
    );
}

// The id of a record named by its accession's molecule: `ref|NM_000002.1|` is a nucleotide,
// whose `E`, `F` and `x` are no residues (fixture `ids_ps_int`, last record).
#[test]
fn the_molecule_of_the_id_decides_the_residues() {
    let input = b">ref|NM_000002.1|\nACGTEFacgtxACGT\n";
    assert_eq!(
        run_words(&with_algo(&["-outfmt", "interval", "-parse_seqids"]), input),
        (0, "4 - 7\n".to_string(), String::new())
    );
    assert_eq!(
        run_words(&with_algo(&["-outfmt", "interval"]), input),
        (0, "6 - 10\n".to_string(), String::new())
    );
}

// NCBI writes `#` for each byte of `-masking_options` outside the visible ASCII range of the
// ASN.1 VisibleString and posts `Error: (803.7) Blast-db-mask-info.algo-options` with the byte's
// line or offset (convert2blastmask 2.17.0, `-masking_options 'é window=64'`, exit status 0);
// LOSAT's serial writer refuses the byte, as in an archive. The XML form writes the bytes
// (fixture `opt_utf8_xml`).
#[test]
fn a_masking_option_outside_the_visible_ascii_range_is_refused_in_asn() {
    for outfmt in ["maskinfo_asn1_text", "maskinfo_asn1_bin"] {
        let (status, out, err) = run_words(
            &[
                "-masking_algorithm",
                "dust",
                "-masking_options",
                "é window=64",
                "-outfmt",
                outfmt,
            ],
            b">a\nacgt\n",
        );
        assert_eq!((status, out), (1, String::new()), "{outfmt}");
        assert_eq!(
            err,
            "Error: an object on which LOSAT's serial writer stops (a character 0xc3 outside the \
             visible ASCII range in an ASN.1 VisibleString is not supported by LOSAT) is not \
             supported by LOSAT's CONVERT2BLASTMASK\n",
            "{outfmt}"
        );
    }
}
