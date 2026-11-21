//! The names outfmt 5 writes: the query and subject definition lines, the `Iteration_message`,
//! and the refusal of what NCBI crashes on (DESIGN_OF4 1.5, 2.2, 2.3). NCBI's XML stream writes
//! every other text itself, a byte >= 0x80 or a control byte included (`report::serial::xml`).

use crate::report::defline::generate_defline;

/// NCBI's `-outfmt 5` crashes making the text; LOSAT refuses the search.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Refusal {
    /// `x_CleanAndCompress` reads past the end of the title (NCBI crashes there).
    ReadsPastEnd,
}

impl Refusal {
    /// The refusal as LOSAT's error: "`<what>` with `<reason>` is not supported by LOSAT's `<PROGRAM>`".
    pub fn error(self, program: &str, what: &str) -> anyhow::Error {
        anyhow::anyhow!(
            "{what} with {} is not supported by LOSAT's {program}",
            self.describe()
        )
    }

    /// The reason, as the end of "a query title with ..." in a message.
    pub fn describe(self) -> &'static str {
        match self {
            Refusal::ReadsPastEnd => {
                "punctuation that NCBI BLAST+ 2.17.0 crashes on when it makes the definition line"
            }
        }
    }
}

/// `CDeflineGenerator().GenerateDefline(bh)` of a record with the FASTA title `title` (no
/// flags): the definition line, or `None` where it throws (`NStr::HtmlDecode` cannot guess the
/// title's encoding, `report::defline::UnknownEncoding`), which NCBI's XML writers catch
/// (`QueryText::id`, `SubjectNames::of_unknown`).
fn generated(title: &[u8], protein: bool) -> Result<Option<Vec<u8>>, Refusal> {
    let Ok(defline) = generate_defline(title, protein, false) else {
        return Ok(None);
    };
    if defline.reads_past_end {
        return Err(Refusal::ReadsPastEnd);
    }
    Ok(Some(defline.text))
}

/// What outfmt 5 writes for an empty query or subject title (blastxml_format.cpp:557-565);
/// XML2 omits the query title and writes an empty subject title.
pub const NO_DEFINITION_LINE: &[u8] = b"No definition line";

/// The definition line of a query as generated (`GenerateDefline`, blastxml_format.cpp:
/// 557-565; blastxml2_format.cpp:590-592): empty for an empty nucleotide title, which
/// `Iteration_query-def` replaces by `No definition line` and `Search_query-title` omits; `None`
/// where `GenerateDefline` throws.
pub fn xml_query_title(title: &[u8], protein: bool) -> Result<Option<Vec<u8>>, Refusal> {
    generated(title, protein)
}

/// `xml_query_title` of a query whose ids `-parse_deflines` read: the parts of the title that
/// come from the Seq-ids (`GenerateDefline`, `report::defline::BioseqFacts`).
fn believed_query_title(
    title: &[u8],
    protein: bool,
    believed: &crate::blastinput::fasta_reader::believed::BelievedIds,
) -> Result<Option<Vec<u8>>, Refusal> {
    let ids = crate::report::db_subject::bioseq_facts(&believed.ids);
    let Ok(defline) = crate::report::defline::generate_bioseq_defline(title, protein, false, &ids)
    else {
        return Ok(None);
    };
    if defline.reads_past_end {
        return Err(Refusal::ReadsPastEnd);
    }
    Ok(Some(defline.text))
}

/// `kSeqId.AsFastaString()` of a query or `-subject` record: `lcl|` and the local id, or the
/// FASTA form of the best id that `-parse_deflines` read; what NCBI's XML writers name a record
/// whose `GenerateDefline` threw.
pub fn fasta_id(
    local_id: &[u8],
    believed: Option<&crate::blastinput::fasta_reader::believed::BelievedIds>,
) -> Vec<u8> {
    let mut out = Vec::new();
    match believed {
        Some(ids) => crate::seqid::label::fasta(&ids.best, &mut out),
        None => {
            out.extend_from_slice(b"lcl|");
            out.extend_from_slice(local_id);
        }
    }
    out
}

/// The `Hit_id` and the title of a `-subject` record: the first word of the generated
/// definition line (`CShowBlastDefline::GetSeqIdList`, showdefline.cpp:218-242), or the local
/// ID `Subject_N` when that is empty; the definition line as generated (empty for an empty
/// nucleotide title).
pub fn xml_subject_names(
    local_id: &[u8],
    title: &[u8],
    protein: bool,
) -> Result<Option<(Vec<u8>, Vec<u8>)>, Refusal> {
    let Some(text) = generated(title, protein)? else {
        return Ok(None);
    };
    let word = text.split(|&byte| byte == b' ').next().unwrap_or(&[]);
    let id = if word.is_empty() {
        local_id.to_vec()
    } else {
        word.to_vec()
    };
    Ok(Some((id, text)))
}

/// The `Iteration_message` of a query (data4xmlformat.cpp:156-169), or `None` when empty.
/// `query_id` is the id that NCBI's warnings carry (`report::query_warnings::warning_query_id`),
/// `warnings` the query's warning messages (any order; sorted and without repeats as NCBI keeps
/// them) and `has_hits` whether the query has Seq-aligns.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:156-169
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:258-293
/// (`GetErrorStrings` is the id and `: ` alone when the query has messages, all warnings in
/// LOSAT.)
pub fn xml_message(query_id: &[u8], warnings: &[String], has_hits: bool) -> Option<Vec<u8>> {
    let mut message: Vec<u8> = Vec::new();
    if !warnings.is_empty() {
        // GetErrorStrings(): the id and ": "
        message.extend_from_slice(query_id);
        message.extend_from_slice(b": ");
        // " " and GetWarningStrings()
        message.push(b' ');
        message.extend_from_slice(query_id);
        message.extend_from_slice(b": ");
        for warning in crate::report::query_warnings::sorted_query_messages(warnings) {
            message.extend_from_slice(warning.as_bytes());
            message.push(b' ');
        }
    }
    if !has_hits {
        if !message.is_empty() {
            message.push(b' ');
        }
        message.extend_from_slice(b"No hits found");
    }
    (!message.is_empty()).then_some(message)
}

/// The title and message of a query, made before the document is written.
pub struct QueryText {
    /// The generated title (`xml_query_title`; empty where `GenerateDefline` threw).
    pub title: Vec<u8>,
    pub message: Option<Vec<u8>>,
    /// Where `GenerateDefline` threw: the query's id that NCBI's catch writes in place of the
    /// record's (`kSeqId.AsFastaString()`, `fasta_id`), with no title (outfmt 5: `No
    /// definition line`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:552-565
    /// (The header, blastxml_format.cpp:700-711, and XML2's `Search`, blastxml2_format.cpp:
    /// 585-595, catch it the same way.)
    pub id: Option<Vec<u8>>,
}

impl QueryText {
    /// The query's id in the document: `record_id`, or the one NCBI's catch writes (`id`).
    pub fn query_id<'a>(&'a self, record_id: &'a [u8]) -> &'a [u8] {
        self.id.as_deref().unwrap_or(record_id)
    }
}

/// The `QueryText` of one query (`program` names the program in a refusal): its record's
/// `local_id`, `warning_id` (`FastaRecord::warning_id`: `Query_N` without `-parse_deflines`) and
/// FASTA `title`, the query's warning messages and whether it has Seq-aligns.
#[allow(clippy::too_many_arguments)]
pub fn query_text(
    program: &str,
    local_id: &[u8],
    warning_id: &[u8],
    believed: Option<&crate::blastinput::fasta_reader::believed::BelievedIds>,
    title: &[u8],
    protein: bool,
    warnings: &[String],
    has_hits: bool,
) -> anyhow::Result<QueryText> {
    let generated = match believed {
        Some(ids) => believed_query_title(title, protein, ids),
        None => xml_query_title(title, protein),
    }
    .map_err(|refusal| refusal.error(program, "a query title"))?;
    let message_id = crate::report::query_warnings::warning_query_id(warning_id, title);
    let message = xml_message(&message_id, warnings, has_hits);
    Ok(match generated {
        Some(title) => QueryText {
            title,
            message,
            id: None,
        },
        None => QueryText {
            title: Vec::new(),
            message,
            id: Some(fasta_id(local_id, believed)),
        },
    })
}

/// Whether NCBI's XML stream throws in the outfmt 5 `Iteration` of a query before its first `Hit`
/// (the query's record as `query_text` takes it): at its `Iteration_query-ID` or
/// `Iteration_query-def` (`report::serial::xml::us_ascii_stream_throws`), the strings it writes
/// before the hits. NCBI then makes none of the loads of the report's reads, whose texts go to its
/// next loads that make the same asks (`report::query_warnings::ReadCursor::cut`). A title that
/// LOSAT refuses (`Refusal`) is not written, so it cuts nothing.
pub fn iteration_cut_before_hits(
    local_id: &[u8],
    warning_id: &[u8],
    believed: Option<&crate::blastinput::fasta_reader::believed::BelievedIds>,
    title: &[u8],
    protein: bool,
) -> bool {
    use crate::report::serial::xml::us_ascii_stream_throws;
    query_text(
        "",
        local_id,
        warning_id,
        believed,
        title,
        protein,
        &[],
        true,
    )
    .is_ok_and(|text| {
        us_ascii_stream_throws(text.query_id(local_id)) || us_ascii_stream_throws(&text.title)
    })
}

/// The `QueryText` of every query of `records` (`program` names the program in a refusal):
/// `warnings[q]` are the query's warning messages and `has_hits(q)` whether it has Seq-aligns.
pub fn query_texts(
    program: &str,
    records: &[crate::blastinput::fasta_reader::FastaRecord],
    protein: bool,
    warnings: &[Vec<String>],
    has_hits: impl Fn(usize) -> bool,
) -> anyhow::Result<Vec<QueryText>> {
    records
        .iter()
        .enumerate()
        .map(|(q_idx, record)| {
            query_text(
                program,
                record.local_id.as_bytes(),
                record.warning_id(),
                record.believed(),
                &record.title,
                protein,
                &warnings[q_idx],
                has_hits(q_idx),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
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
        let o = "One or more O characters replaced by X for alignment score calculations at positions 20";
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
}
