//! Warnings that NCBI writes to standard error for a search, whatever the output format.

use bio::io::fasta;
use std::io::{self, Write};

/// The warnings that NCBI posts between the reports of a search: the title warnings of a
/// query batch when it is read (before the report of its first query) and the warning of
/// an invalid query when its report starts (`PrintOneResultSet`, before the query's
/// preamble). NCBI posts them on `cerr`, which the C++ library ties to `cout`, so the
/// report written so far reaches standard output before each of them; a failed write of
/// the report stops NCBI there (its stream throws, `blast_format.cpp:118-119`).
///
/// NCBI reference: ncbi-blast/c++/src/corelib/ncbidiag.cpp:4086
/// ```c
///         CDiagHandler* handler = new CStreamDiagHandler(&NcbiCerr, true, kLogName_Stderr);
/// ```
/// NCBI reference: ncbi-blast/c++/include/corelib/ncbistre.hpp:543-544
/// ```c
/// #define NcbiCout                 IO_PREFIX::cout
/// #define NcbiCerr                 IO_PREFIX::cerr
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1450-1452
/// ```c
///     if (results.HasWarnings()) {
///         ERR_POST(Warning << results.GetWarningStrings());
///     }
/// ```
pub struct QueryWarnings<'a> {
    /// The warnings to write before the report of each query, by query index.
    pub before: &'a [Vec<u8>],
    /// Standard error (the CLI) or the caller's diagnostics.
    pub sink: &'a mut (dyn Write + Send),
}

impl QueryWarnings<'_> {
    /// Writes the warnings of query `index`, after the report written so far (`report`)
    /// is flushed; a failed flush returns before them.
    pub fn before_query(&mut self, index: usize, report: &mut dyn Write) -> io::Result<()> {
        match self.before.get(index) {
            Some(warnings) if !warnings.is_empty() => {
                report.flush()?;
                self.sink.write_all(warnings)?;
                self.sink.flush()
            }
            _ => Ok(()),
        }
    }
}

/// The warning for a query whose ungapped Karlin-Altschul parameters cannot be computed.
///
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:2783-2790
/// ```c
///       if (loop_status) {
///           contexts[context].is_valid = FALSE;
///           ...
///           if (!Blast_QueryIsTranslated(program) ) {
///              Blast_MessageWrite(blast_message, eBlastSevWarning, context,
///              kBlastErrMsg_CantCalculateUngappedKAParams);
///           }
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_message.c:37-40
/// ```c
/// const char* kBlastErrMsg_CantCalculateUngappedKAParams
///     = "Could not calculate ungapped Karlin-Altschul parameters due "
///       "to an invalid query sequence or its translation. Please verify the "
///       "query sequence(s) and/or filtering options";
/// ```
/// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:534-543
/// ```c++
///                 const string kTitle = queries.GetTitle(index);
///                 string query_id = id->GetSeqIdString();
///                 if (kTitle != kEmptyStr) {
///                     query_id += " " + kTitle;
///                 }
///                  if(query_id.size() > 35) {
///                 	 query_id = query_id.substr(0, 25) + ".. ";
///                  }
///
///                 messages[index].SetQueryId(query_id);
/// ```
/// NCBI reference: c++/src/algo/blast/api/blast_results.cpp:283-291
/// ```c++
///     string retval(m_Errors.GetQueryId());
///     if ( !retval.empty() ) {    // in case the query id is not known
///         retval += ": ";
///     }
///     ITERATE(TQueryMessages, iter, m_Errors) {
///         if ((**iter).GetSeverity() == eBlastSevWarning) {
///             retval += (*iter)->GetMessage(false) + " ";
///         }
///     }
/// ```
/// The local ID of an input query is `Query_<n>` and its title is the FASTA defline.
///
/// `index` is the 0-based position of the query in the input; `program` is the lower-case
/// program name that NCBI's diagnostics print (`[tblastn]`, `[blastn]`).
pub fn invalid_query_warning(program: &str, index: usize, query: &fasta::Record) -> Vec<u8> {
    query_warning(program, index, query, &[INVALID_QUERY_MESSAGE.to_string()])
}

/// NCBI's `kBlastErrMsg_CantCalculateUngappedKAParams` (see `invalid_query_warning`).
pub const INVALID_QUERY_MESSAGE: &str = "Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options";

/// The warning line of query `index` with its warning messages, each followed by a space
/// (`GetWarningStrings`, see `invalid_query_warning`); empty when there is none.
pub fn query_warning(
    program: &str,
    index: usize,
    query: &fasta::Record,
    messages: &[String],
) -> Vec<u8> {
    if messages.is_empty() {
        return Vec::new();
    }
    let mut query_id = format!("Query_{} {}", index + 1, query.id()).into_bytes();
    if let Some(desc) = query.desc() {
        query_id.extend_from_slice(b" ");
        query_id.extend_from_slice(desc.as_bytes());
    }
    if query_id.len() > 35 {
        query_id.truncate(25);
        query_id.extend_from_slice(b".. ");
    }
    let mut warning = format!("Warning: [{program}] ").into_bytes();
    warning.extend_from_slice(&query_id);
    warning.extend_from_slice(b": ");
    for message in messages {
        warning.extend_from_slice(message.as_bytes());
        warning.push(b' ');
    }
    warning.push(b'\n');
    warning
}

/// The warning of a protein query with pyrrolysine (O), which NCBI reads as X, or `None`.
///
/// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:920-932
/// ```c
///     if (warnings && replaced_residues.size() > 0) {
///         *warnings += "One or more O characters replaced by X for ";
///         *warnings += "alignment score calculations at positions ";
///         *warnings += NStr::IntToString(replaced_residues[0]);
///         for (i = 1; i < min(kMaxResiduesToWarnAbout, replaced_residues.size());
///              i++) {
///             *warnings += ", " + NStr::IntToString(replaced_residues[i]);
///         }
///         if (replaced_residues.size() > kMaxResiduesToWarnAbout) {
///             *warnings += ",... (only first ";
///             *warnings += NStr::SizetToString(kMaxResiduesToWarnAbout);
///             *warnings += " shown)";
///         }
///     }
/// ```
/// The query's O messages come before its Karlin-Altschul message (blast_setup_cxx.cpp:
/// 608-621 adds them when the query is read).
pub fn replaced_o_message(sequence: &[u8]) -> Option<String> {
    let positions: Vec<usize> = sequence
        .iter()
        .enumerate()
        .filter(|(_, residue)| residue.eq_ignore_ascii_case(&b'O'))
        .map(|(position, _)| position)
        .collect();
    let first = *positions.first()?;
    let mut message = format!(
        "One or more O characters replaced by X for alignment score calculations at positions {first}"
    );
    for position in positions.iter().take(20).skip(1) {
        message.push_str(&format!(", {position}"));
    }
    if positions.len() > 20 {
        message.push_str(",... (only first 20 shown)");
    }
    Some(message)
}

/// The warning for a hit list size below 5.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
/// ```c++
///     if(hitlist_size < 5){
///    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
///     }
/// ```
pub fn few_matches_warning(program: &str) -> Vec<u8> {
    format!("Warning: [{program}] Examining 5 or more matches is recommended\n").into_bytes()
}

#[cfg(test)]
mod tests {
    use super::*;

    // Query_1 + FASTA title is shortened after 35 bytes, then the invalid-Karlin
    // warning keeps NCBI's trailing space and newline.
    #[test]
    fn invalid_query_warning_matches_ncbi_bytes() {
        let short = fasta::Record::with_attrs("nohit_query", None, b"W");
        assert_eq!(
            invalid_query_warning("tblastn", 0, &short),
            b"Warning: [tblastn] Query_1 nohit_query: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options \n"
        );
        let long =
            fasta::Record::with_attrs("long_header", Some("abcdefghijklmnopqrstuvwxyz"), b"W");
        assert!(invalid_query_warning("blastn", 1, &long)
            .starts_with(b"Warning: [blastn] Query_2 long_header abcde.. : "));
    }
}
