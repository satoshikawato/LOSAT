use super::*;

const KA: &str = "Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options";

#[test]
fn messages_of_the_probes() {
    // blastp_sub_invalid: two spaces at both joints
    let message = xml_message(b"Query_1 inv all X", &[KA.to_string()], false).unwrap();
    assert_eq!(
        String::from_utf8(message).unwrap(),
        format!("Query_1 inv all X:  Query_1 inv all X: {KA}  No hits found")
    );
    // with hits: the trailing space of the warning ends the message
    let o =
        "One or more O characters replaced by X for alignment score calculations at positions 20";
    let message = xml_message(b"Query_1 qp1 rare letters", &[o.to_string()], true).unwrap();
    assert_eq!(
        String::from_utf8(message).unwrap(),
        format!("Query_1 qp1 rare letters:  Query_1 qp1 rare letters: {o} ")
    );
    // no message and hits: no element; no hits alone: "No hits found" without a space
    assert_eq!(xml_message(b"Query_1", &[], true), None);
    assert_eq!(
        xml_message(b"Query_1", &[], false).unwrap(),
        b"No hits found"
    );
    // messages are sorted and without repeats
    let message = xml_message(b"Q", &["b".into(), "a".into(), "b".into()], true).unwrap();
    assert_eq!(String::from_utf8(message).unwrap(), "Q:  Q: a b ");
}

#[test]
fn an_iteration_is_cut_before_its_hits_by_a_title_outside_us_ascii() {
    // K4h review 3 H1 (E2): q1's two-byte UTF-8 title cuts its Iteration at the query-def,
    // before any Hit; an ASCII title cuts nothing.
    assert!(iteration_cut_before_hits(
        b"Query_2",
        b"Query_2",
        None,
        b"q1 a\xc3\xa9b z",
        true
    ));
    assert!(!iteration_cut_before_hits(
        b"Query_2", b"Query_2", None, b"q1 ab z", true
    ));
    // A title whose GenerateDefline throws (0x81, in no encoding NCBI guesses) is not written:
    // the catch's id `lcl|Query_2` and no title come before the hits.
    assert!(!iteration_cut_before_hits(
        b"Query_2",
        b"Query_2",
        None,
        b"q1 a\x81b",
        true
    ));
}

#[test]
fn query_and_subject_names_of_the_probe_titles() {
    // blastp_sub_titles (DESIGN_OF4 2.2, 2.3)
    assert_eq!(
        xml_query_title(b"q1 title ending with periods...", true).unwrap(),
        Some(b"q1 title ending with periods".to_vec())
    );
    assert_eq!(
        xml_query_title(b"TPA: q3 prefixed", true).unwrap(),
        Some(b"q3 prefixed".to_vec())
    );
    assert_eq!(
        xml_query_title(b"q4  double  space ,  comma ;", true).unwrap(),
        Some(b"q4 double space, comma".to_vec())
    );
    assert_eq!(
        xml_query_title(b"", true).unwrap(),
        Some(b"unnamed protein product".to_vec())
    );
    // an empty nucleotide title stays empty: outfmt 5 writes `No definition line` for it,
    // XML2 omits the member
    assert_eq!(xml_query_title(b"", false).unwrap(), Some(Vec::new()));
    let (id, def) = xml_subject_names(b"Subject_1", b"s1 subj ending.", true)
        .unwrap()
        .unwrap();
    assert_eq!(
        (id.as_slice(), def.as_slice()),
        (&b"s1"[..], &b"s1 subj ending"[..])
    );
    let (id, def) = xml_subject_names(b"Subject_2", b"s2 &amp; <y> \"q\"", true)
        .unwrap()
        .unwrap();
    assert_eq!(id, b"s2");
    // decoded; the encoder escapes it again
    assert_eq!(def, b"s2 & <y> \"q\"");
    let (id, def) = xml_subject_names(b"Subject_4", b"TPA: s4 tpa", true)
        .unwrap()
        .unwrap();
    assert_eq!(
        (id.as_slice(), def.as_slice()),
        (&b"s4"[..], &b"s4 tpa"[..])
    );
    // an empty protein title: `unnamed`; an empty nucleotide title: the local ID
    let (id, def) = xml_subject_names(b"Subject_5", b"", true).unwrap().unwrap();
    assert_eq!(
        (id.as_slice(), def.as_slice()),
        (&b"unnamed"[..], &b"unnamed protein product"[..])
    );
    let (id, def) = xml_subject_names(b"Subject_1", b"", false)
        .unwrap()
        .unwrap();
    assert_eq!(
        (id.as_slice(), def.as_slice()),
        (&b"Subject_1"[..], &b""[..])
    );
}

#[test]
fn titles_ncbi_cannot_decode_are_its_catch_and_the_rest_reaches_the_stream() {
    // K4h probe q.mq classes: valid UTF-8 and Latin-1 titles are written by the stream
    // (which throws there), an undecidable encoding is NCBI's catch (`None`). A last byte
    // >= 0x80 is dropped by `x_CleanAndCompress` (u2end: `q33 one\xc3`).
    assert_eq!(
        xml_query_title("caf\u{e9} x".as_bytes(), true).unwrap(),
        Some("caf\u{e9} x".as_bytes().to_vec())
    );
    assert_eq!(
        xml_query_title("caf\u{e9}".as_bytes(), true).unwrap(),
        Some(b"caf\xc3".to_vec())
    );
    assert_eq!(
        xml_query_title(b"q6 a\xe9b z", true).unwrap(),
        Some(b"q6 a\xc3\xa9b z".to_vec())
    );
    assert_eq!(xml_query_title(b"q9 a\x81b z", true).unwrap(), None);
    assert_eq!(
        xml_subject_names(b"Subject_8", b"s8 a\x81b z", true),
        Ok(None)
    );
    assert_eq!(fasta_id(b"Query_9", None), b"lcl|Query_9");
    let text = query_text(
        "BLASTP",
        b"Query_9",
        b"Query_9",
        None,
        b"q9 a\x81b z",
        true,
        &[],
        true,
    )
    .unwrap();
    assert_eq!(
        (text.title.as_slice(), text.query_id(b"Query_9")),
        (&b""[..], &b"lcl|Query_9"[..])
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-719
/// A result shown with an earlier query's id, title and length (`labelled_query_text`, line RG-KB
/// item 7c) has its own message and no title: NCBI's XML formatter generates no definition line
/// of its record, so a title that `query_text` refuses (`Refusal::ReadsPastEnd`) is not refused.
#[test]
fn a_labelled_query_has_its_message_and_no_title() {
    let refused = [&b","[..], b", ,", b" ,", b";", b"a ,"]
        .into_iter()
        .find(|title| xml_query_title(title, true).is_err())
        .expect("a title that NCBI crashes on");
    let own = query_text(
        "BLASTP",
        b"Query_2",
        b"Query_2",
        None,
        refused,
        true,
        &[],
        false,
    );
    assert!(own.is_err());
    let text = labelled_query_text(b"Query_2", refused, &[], false);
    assert!(text.title.is_empty() && text.id.is_none());
    assert_eq!(text.message.as_deref(), Some(&b"No hits found"[..]));
}
