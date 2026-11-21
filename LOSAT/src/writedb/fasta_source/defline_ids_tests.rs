use super::*;

#[test]
fn ranges_of_the_id_token() {
    assert_eq!(parse_range(b"seqA:1-10"), (5, 0, 9));
    assert_eq!(parse_range(b"seqB:c10-1"), (6, 9, 0));
    assert_eq!(parse_range(b":1-10"), (5, 0, 9));
    assert_eq!(parse_range(b"x:10-1"), (0, 0, 0));
    assert_eq!(parse_range(b"x:c1-10"), (0, 0, 0));
    assert_eq!(parse_range(b"x:1-"), (0, 0, 0));
    assert_eq!(parse_range(b"x:-1"), (0, 0, 0));
    assert_eq!(parse_range(b"1-10"), (0, 0, 0));
    assert_eq!(parse_range(b"NC_000914.1:1-50"), (5, 0, 49));
    assert_eq!(parse_range(b"a:5-5"), (4, 4, 4));
    assert_eq!(parse_range(b""), (0, 0, 0));
}

fn check_with(
    id: &str,
    protein: bool,
    quick_id_check: bool,
) -> (Result<Vec<ParsedSeqId>, ParseError>, Vec<String>) {
    let mut warnings = Vec::new();
    let result = process_ids(
        id.as_bytes(),
        protein,
        quick_id_check,
        true,
        7,
        &mut |w| match w {
            ReaderWarning::Line(line) => warnings.push(String::from_utf8_lossy(line).into_owned()),
            ReaderWarning::SeqId(w) => warnings.push(format!(
                "(1306.{}) {}",
                w.subcode,
                String::from_utf8_lossy(&w.message)
            )),
        },
    );
    (result, warnings)
}

fn check(id: &str, protein: bool) -> (Result<Vec<ParsedSeqId>, String>, Vec<String>) {
    let (result, warnings) = check_with(id, protein, true);
    (result.map_err(|e| e.message), warnings)
}

// Messages of makeblastdb 2.17.0 (survey B fixtures p05_raw, p06*, p07, p08, n03, n04).
#[test]
fn checks_and_messages_of_the_first_parse() {
    let (ids, warnings) = check("a,b", true);
    assert_eq!(
        ids.unwrap(),
        [ParsedSeqId::Local(ObjectId::Str(b"a_b".to_vec()))]
    );
    assert_eq!(
        warnings,
        ["Near line 7, the sequence id string contains 'comma' symbol, which has been replaced with 'underscore' symbol. Please correct the sequence id string."]
    );
    let acc31 = format!("ref|NP_{}|", "1".repeat(28));
    assert_eq!(
        check(&acc31, true).0.unwrap_err(),
        "Near line 7, the accession is too long.  Its length is 31 but the maximum allowed accession length is 30.  Please find and correct all accessions that are too long."
    );
    assert!(check(&format!("ref|NP_{}|", "1".repeat(27)), true)
        .0
        .is_ok());
    let local51 = format!("lcl|{}", "a".repeat(51));
    let (result, warnings) = check(&local51, true);
    assert_eq!(
        result.unwrap_err(),
        "Near line 7, the local id is too long.  Its length is 51 but the maximum allowed local id length is 50.  Please find and correct all local ids that are too long."
    );
    assert_eq!(
        warnings,
        ["Fasta Reader: sequence id ends with 51 valid amino-acid characters.  Was the sequence accidentally placed in the definition line?"]
    );
    assert_eq!(
        check("lcl|[abc]", true).0.unwrap_err(),
        "'[abc]' is not a valid local ID"
    );
    assert!(check("lcl|a[bc]", true).0.is_ok());
    let nuc21 = format!("lcl|x{}", "ACGT".repeat(5) + "A");
    let (result, warnings) = check(&nuc21, false);
    assert!(result.is_ok());
    assert_eq!(
        warnings,
        ["Fasta Reader: sequence id ends with 21 valid nucleotide characters.  Was the sequence accidentally placed in the definition line?"]
    );
    let nuc26 = format!("lcl|x{}", "acgt".repeat(6) + "AC");
    assert_eq!(
        check(&nuc26, false).0.unwrap_err(),
        "Fasta Reader: sequence id ends with 26 valid nucleotide characters.  Was the sequence accidentally placed in the definition line?"
    );
    assert_eq!(
        check("gi|abc|", true).0.unwrap_err(),
        "Could not construct seq-id from 'gi|abc|'"
    );
}

// dustmasker 2.17.0 -parse_seqids (`tests/dustmasker/reader`, cases ps_badlocal*,
// ps_nucid_err, ps_seqidwarn) and segmasker (seg_longid_ps): without `fQuickIDCheck`
// every byte of a local id is checked; the codes and positions of the exceptions.
#[test]
fn the_full_id_check_and_the_codes_of_the_errors() {
    assert!(check_with("a]b", false, true).0.is_ok());
    let error = check_with("a]b", false, false).0.unwrap_err();
    assert_eq!(error.message, "'a]b' is not a valid local ID");
    assert_eq!(
        (error.code, error.pos, error.thrower),
        (ErrCode::InvalidId, 7, Thrower::PostError)
    );
    assert!(check_with("x[y]", true, false).0.is_err());
    assert!(check_with("lcl|a\"b", true, false).0.is_err());
    let error = check_with(&"a".repeat(51), true, false).0.unwrap_err();
    assert_eq!(error.code, ErrCode::IdTooLong);
    let error = check_with(&"a".repeat(51), false, false).0.unwrap_err();
    assert_eq!(error.code, ErrCode::Format);
    let (result, warnings) = check_with("xyz|ABC|DEF", false, false);
    let error = result.unwrap_err();
    assert_eq!(
        error.message,
        "Could not construct seq-id from 'xyz|ABC|DEF'"
    );
    assert_eq!(error.code, ErrCode::NoIds);
    assert_eq!(warnings.len(), 2, "{warnings:?}");
}
