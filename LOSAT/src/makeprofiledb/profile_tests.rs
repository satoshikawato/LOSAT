use super::*;

/// The fixture `fc.smp` (frequency ratios only) with its query's `seq-data` replaced.
fn profile(seq_data: &str) -> PssmWithParameters {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/makeprofiledb/in/fc.smp"
    ))
    .unwrap();
    let start = text.find("seq-data ncbistdaa '").unwrap();
    let end = start + text[start..].find("'H").unwrap() + 2;
    let text = format!("{}{seq_data}{}", &text[..start], &text[end..]);
    crate::blastinput::pssm_text::read_pssm_with_parameters(text.as_bytes()).unwrap()
}

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1674-1707,1748-1781
// Every byte of the two tables: NCBIeaa's capitals, `*` (25) and `-` (0), IUPACaa's letters of
// either case (the capital's code), 255 for any other byte.
#[test]
fn every_byte_has_ncbis_code() {
    const LETTERS: &[u8] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
    let code_of = |letter: u8| LETTERS.iter().position(|&l| l == letter).map(|c| c as u8);
    let bytes: Vec<u8> = (0..=255).collect();
    let eaa = crate::writedb::convert::ncbieaa_codes(&bytes);
    let iupac = crate::writedb::convert::iupacaa_codes(&bytes);
    for byte in 0..=255u8 {
        let want_eaa = match byte {
            b'A'..=b'Z' | b'*' | b'-' => code_of(byte).unwrap(),
            _ => 255,
        };
        let want_iupac = if byte.is_ascii_alphabetic() {
            code_of(byte.to_ascii_uppercase()).unwrap()
        } else {
            255
        };
        assert_eq!(eaa[usize::from(byte)], want_eaa, "ncbieaa 0x{byte:02x}");
        assert_eq!(iupac[usize::from(byte)], want_iupac, "iupacaa 0x{byte:02x}");
    }
}

// NCBI makeprofiledb 2.17.0 (row PS10c probes, 2026-11-23): IUPACaa of either case builds the
// database of its capitals (fixture row QC01); NCBIeaa other than the capitals, `*` and `-`, and
// IUPACaa other than letters, have the code 255.
#[test]
fn query_codes_follow_ncbis_conversion() {
    let iupacaa = profile("seq-data iupacaa \"kDwEkl\"");
    assert_eq!(query_ncbistdaa(&iupacaa).unwrap(), [10, 4, 20, 5, 10, 11]);
    let iupacaa = profile("seq-data iupacaa \"*-K\"");
    assert_eq!(query_ncbistdaa(&iupacaa).unwrap(), [255, 255, 10]);
    let ncbieaa = profile("seq-data ncbieaa \"*-+ka@[K\"");
    assert_eq!(
        query_ncbistdaa(&ncbieaa).unwrap(),
        [25, 0, 255, 255, 255, 255, 255, 10]
    );
    let ncbistdaa = profile("seq-data ncbistdaa '1EFF0A'H");
    assert_eq!(query_ncbistdaa(&ncbistdaa).unwrap(), [30, 255, 10]);
}

// NCBI ends with SIGSEGV on such a query of the right length (fixture rows QC02-QC05: a length
// that does not match, or a matrix the engine does not know, is NCBI's error first).
#[test]
fn a_code_past_the_residues_is_refused_where_the_engine_counts_it() {
    let refusal = |data: &str, matrix: &str| {
        let pssm = profile(data);
        let query = query_ncbistdaa(&pssm).unwrap();
        query_code_past_residues(&pssm, &query, matrix).map(|error| match error {
            ProfileError::Rps(RpsError::Refused(text)) => text,
            other => panic!("{other:?}"),
        })
    };
    let text = refusal("seq-data ncbieaa \"KD+E\"", "BLOSUM62").unwrap();
    assert!(text.contains("the letter 0x2b at position 3"), "{text}");
    let text = refusal("seq-data iupacaa \"kd*\"", "PAM30").unwrap();
    assert!(text.contains("the letter 0x2a at position 3"), "{text}");
    for code in [28, 29, 30, 31, 32, 255] {
        let text = refusal(&format!("seq-data ncbistdaa '0A{code:02X}'H"), "BLOSUM62").unwrap();
        assert!(
            text.contains(&format!("the NCBIstdaa code {code} at position 2")),
            "{text}"
        );
    }
    assert_eq!(refusal("seq-data ncbieaa \"KD+E\"", "FOO"), None);
    assert_eq!(refusal("seq-data iupacaa \"kd\"", "BLOSUM62"), None);
    assert_eq!(refusal("seq-data ncbistdaa '0A1B'H", "BLOSUM62"), None);
}

// NCBI makeprofiledb 2.17.0 (row PS10d survey §2, 2026-11-23): the query of a profile with
// frequency ratios only in a nucleotide coding stops the build with `CSeqConvert`'s exception
// (`Error: ERROR: Unknown conversion.`, exit status 255) in each coding; empty data is an empty
// query, whose length NCBI then finds wrong.
#[test]
fn a_nucleotide_query_is_ncbis_unknown_conversion() {
    for data in [
        "seq-data iupacna \"ACGTACGTACGT\"",
        "seq-data ncbi2na '1B1B1B'H",
        "seq-data ncbi4na '124812481248'H",
    ] {
        match query_ncbistdaa(&profile(data)) {
            Err(ProfileError::Exception(exception)) => {
                assert_eq!(exception.get_msg(), b"Unknown conversion.", "{data}")
            }
            other => panic!("{data}: {other:?}"),
        }
    }
    for data in [
        "seq-data iupacna \"\"",
        "seq-data ncbi2na ''H",
        "seq-data ncbi4na ''H",
    ] {
        assert_eq!(query_ncbistdaa(&profile(data)).unwrap(), [0u8; 0], "{data}");
    }
}
