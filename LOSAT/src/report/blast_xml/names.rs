//! The names outfmt 5 writes: the query and subject definition lines, the `Iteration_message`,
//! and the refusal of what NCBI cannot write (DESIGN_OF4 1.5, 2.2, 2.3).

use crate::report::defline::generate_defline;

/// NCBI's `-outfmt 5` cannot write the text (or crashes making it); LOSAT refuses the search.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Refusal {
    /// A byte >= 0x80: the US-ASCII stream throws (`CUtf8::SymbolToChar`, corelib/ncbistr.cpp:
    /// 6900-6903): rc 255 for the first query's header, a failed Iteration or Hit otherwise.
    NonAscii,
    /// A control byte (0x01-0x1F other than TAB, LF, CR), which NCBI's stream replaces by `#`
    /// with an error message.
    Control,
    /// `GenerateDefline` throws: `NStr::HtmlDecode` cannot guess the title's encoding
    /// (the title has a byte >= 0x80).
    UnknownEncoding,
    /// `x_CleanAndCompress` reads past the end of the title (NCBI crashes there).
    ReadsPastEnd,
}

impl Refusal {
    /// The refusal as LOSAT's error: "<what> with <reason> is not supported by LOSAT's <PROGRAM>".
    pub fn error(self, program: &str, what: &str) -> anyhow::Error {
        anyhow::anyhow!(
            "{what} with {} is not supported by LOSAT's {program}",
            self.describe()
        )
    }

    /// The reason, as the end of "a query title with ..." in a message.
    pub fn describe(self) -> &'static str {
        match self {
            Refusal::NonAscii => "a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write)",
            Refusal::Control => "a control character (which NCBI BLAST+ replaces by '#' with an error message in the XML output)",
            Refusal::UnknownEncoding => "an encoding that NCBI BLAST+ cannot guess",
            Refusal::ReadsPastEnd => "punctuation that NCBI BLAST+ 2.17.0 crashes on when it makes the definition line",
        }
    }
}

/// Whether `text` can be written by NCBI's XML stream: no byte >= 0x80, and no control byte
/// other than TAB, LF and CR.
pub fn check_text(text: &[u8]) -> Result<(), Refusal> {
    for &byte in text {
        if byte >= 0x80 {
            return Err(Refusal::NonAscii);
        }
        if byte < 0x20 && !matches!(byte, 0x09 | 0x0a | 0x0d) {
            return Err(Refusal::Control);
        }
    }
    Ok(())
}

/// Refuses a `-db` name that NCBI's XML stream cannot write (a byte outside US-ASCII); each
/// program checks it before its search.
pub fn check_database_name(program: &str, name: &str) -> anyhow::Result<()> {
    check_text(name.as_bytes()).map_err(|refusal| refusal.error(program, "a database name"))
}

/// `CDeflineGenerator().GenerateDefline(bh)` of a record with the FASTA title `title` (no
/// flags), checked for what the XML stream can write.
fn generated(title: &[u8], protein: bool) -> Result<Vec<u8>, Refusal> {
    let defline = generate_defline(title, protein, false).map_err(|_| Refusal::UnknownEncoding)?;
    if defline.reads_past_end {
        return Err(Refusal::ReadsPastEnd);
    }
    check_text(&defline.text)?;
    Ok(defline.text)
}

const NO_DEFINITION_LINE: &[u8] = b"No definition line";

/// `Iteration_query-def` / `BlastOutput_query-def`: the generated definition line, or
/// `No definition line` (blastxml_format.cpp:557-565).
pub fn xml_query_def(title: &[u8], protein: bool) -> Result<Vec<u8>, Refusal> {
    let text = generated(title, protein)?;
    Ok(if text.is_empty() {
        NO_DEFINITION_LINE.to_vec()
    } else {
        text
    })
}

/// The `Hit_id` and `Hit_def` of a `-subject` record: the first word of the generated
/// definition line (`CShowBlastDefline::GetSeqIdList`, showdefline.cpp:218-242), or the local
/// ID `Subject_N` when that is empty; the definition line, or `No definition line`.
pub fn xml_subject_names(
    local_id: &[u8],
    title: &[u8],
    protein: bool,
) -> Result<(Vec<u8>, Vec<u8>), Refusal> {
    let text = generated(title, protein)?;
    let word = text.split(|&byte| byte == b' ').next().unwrap_or(&[]);
    let id = if word.is_empty() {
        local_id.to_vec()
    } else {
        word.to_vec()
    };
    let def = if text.is_empty() {
        NO_DEFINITION_LINE.to_vec()
    } else {
        text
    };
    Ok((id, def))
}

/// The `Hit_def` of a `-db` sequence: its stored title as stored, or `No definition line`
/// (showdefline.cpp:265-303 reads the Blast-def-line's title).
pub fn xml_db_def(stored_title: &[u8]) -> Result<Vec<u8>, Refusal> {
    check_text(stored_title)?;
    Ok(if stored_title.is_empty() {
        NO_DEFINITION_LINE.to_vec()
    } else {
        stored_title.to_vec()
    })
}

/// The `Iteration_message` of a query (data4xmlformat.cpp:156-169), or `None` when empty.
/// `query_id` is the id that NCBI's warnings carry (`report::query_warnings::warning_query_id`),
/// `warnings` the query's warning messages (any order; sorted and without repeats as NCBI keeps
/// them) and `has_hits` whether the query has Seq-aligns.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:156-169
/// ```c++
///                 string errors = results[i].GetErrorStrings();
///                 if (results[i].HasWarnings()) {
///                     if ( !errors.empty() ) {
///                         errors += " ";
///                     }
///                     errors += results[i].GetWarningStrings();
///                 }
///                 if ( !results[i].HasAlignments() ) {
///                     errors += (errors.empty() ? kEmptyStr : " ");
///                     errors += CBlastFormatUtil::kNoHitsFound;
///                 }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:258-293
/// ```c++
///     string retval(m_Errors.GetQueryId());
///     if ( !retval.empty() ) {    // in case the query id is not known
///         retval += ": ";
///     }
///     ITERATE(TQueryMessages, iter, m_Errors) {
///         if ((**iter).GetSeverity() == eBlastSevWarning) {
///             retval += (*iter)->GetMessage(false) + " ";
/// ```
/// (`GetErrorStrings` is the id and `: ` alone when the query has messages, all warnings in
/// LOSAT.)
pub fn xml_message(query_id: &[u8], warnings: &[String], has_hits: bool) -> Option<Vec<u8>> {
    let mut message: Vec<u8> = Vec::new();
    if !warnings.is_empty() {
        let mut sorted: Vec<&String> = warnings.iter().collect();
        sorted.sort();
        sorted.dedup();
        // GetErrorStrings(): the id and ": "
        message.extend_from_slice(query_id);
        message.extend_from_slice(b": ");
        // " " and GetWarningStrings()
        message.push(b' ');
        message.extend_from_slice(query_id);
        message.extend_from_slice(b": ");
        for warning in sorted {
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

/// The definition line and message of a query, made and checked before the document is
/// written.
pub struct QueryText {
    pub def: Vec<u8>,
    pub message: Option<Vec<u8>>,
}

/// The `QueryText` of one query (`program` names the program in a refusal): its record's
/// `local_id` (`Query_N`) and FASTA `title`, the query's warning messages and whether it has
/// Seq-aligns.
pub fn query_text(
    program: &str,
    local_id: &str,
    title: &[u8],
    protein: bool,
    warnings: &[String],
    has_hits: bool,
) -> anyhow::Result<QueryText> {
    let def =
        xml_query_def(title, protein).map_err(|refusal| refusal.error(program, "a query title"))?;
    let id = crate::report::query_warnings::warning_query_id(local_id.as_bytes(), title);
    let message = xml_message(&id, warnings, has_hits);
    if let Some(message) = &message {
        check_text(message).map_err(|refusal| refusal.error(program, "a query title"))?;
    }
    Ok(QueryText { def, message })
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
                &record.local_id,
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
    fn query_and_subject_names_of_the_probe_titles() {
        // blastp_sub_titles (DESIGN_OF4 2.2, 2.3)
        assert_eq!(
            xml_query_def(b"q1 title ending with periods...", true).unwrap(),
            b"q1 title ending with periods"
        );
        assert_eq!(
            xml_query_def(b"TPA: q3 prefixed", true).unwrap(),
            b"q3 prefixed"
        );
        assert_eq!(
            xml_query_def(b"q4  double  space ,  comma ;", true).unwrap(),
            b"q4 double space, comma"
        );
        assert_eq!(
            xml_query_def(b"", true).unwrap(),
            b"unnamed protein product"
        );
        assert_eq!(xml_query_def(b"", false).unwrap(), b"No definition line");
        let (id, def) = xml_subject_names(b"Subject_1", b"s1 subj ending.", true).unwrap();
        assert_eq!(
            (id.as_slice(), def.as_slice()),
            (&b"s1"[..], &b"s1 subj ending"[..])
        );
        let (id, def) = xml_subject_names(b"Subject_2", b"s2 &amp; <y> \"q\"", true).unwrap();
        assert_eq!(id, b"s2");
        // decoded; the encoder escapes it again
        assert_eq!(def, b"s2 & <y> \"q\"");
        let (id, def) = xml_subject_names(b"Subject_4", b"TPA: s4 tpa", true).unwrap();
        assert_eq!(
            (id.as_slice(), def.as_slice()),
            (&b"s4"[..], &b"s4 tpa"[..])
        );
        // an empty protein title: `unnamed`; an empty nucleotide title: the local ID
        let (id, def) = xml_subject_names(b"Subject_5", b"", true).unwrap();
        assert_eq!(
            (id.as_slice(), def.as_slice()),
            (&b"unnamed"[..], &b"unnamed protein product"[..])
        );
        let (id, def) = xml_subject_names(b"Subject_1", b"", false).unwrap();
        assert_eq!(
            (id.as_slice(), def.as_slice()),
            (&b"Subject_1"[..], &b"No definition line"[..])
        );
        // -db: the stored title as stored
        assert_eq!(xml_db_def(b"s1 subj ending.").unwrap(), b"s1 subj ending.");
        assert_eq!(xml_db_def(b"TPA: s4 tpa").unwrap(), b"TPA: s4 tpa");
        assert_eq!(xml_db_def(b"").unwrap(), b"No definition line");
    }

    #[test]
    fn what_ncbi_cannot_write_is_refused() {
        assert!(xml_query_def("caf\u{e9}".as_bytes(), true).is_err());
        assert!(xml_query_def(b"caf\xe9", true).is_err());
        assert_eq!(xml_query_def(b"a\x01b", true), Err(Refusal::Control));
        assert!(xml_subject_names(b"Subject_1", "\u{e9}".as_bytes(), true).is_err());
        assert_eq!(xml_db_def("\u{e9}".as_bytes()), Err(Refusal::NonAscii));
        assert_eq!(check_text(b"tab\there\n"), Ok(()));
    }
}
