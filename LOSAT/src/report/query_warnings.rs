//! Warnings that NCBI writes to standard error for a search, whatever the output format.

use bio::io::fasta;

/// The warning for a query whose ungapped Karlin-Altschul parameters cannot be computed.
///
/// NCBI c++/src/algo/blast/core/blast_stat.c:2780-2792:
/// if (loop_status && !Blast_QueryIsTranslated(program))
///     Blast_MessageWrite(..., eBlastSevWarning, context,
///                        kBlastErrMsg_CantCalculateUngappedKAParams);
/// NCBI c++/src/algo/blast/core/blast_message.c:37-40:
/// kBlastErrMsg_CantCalculateUngappedKAParams = "Could not calculate ...".
/// NCBI c++/src/algo/blast/api/blast_setup_cxx.cpp:535-543:
/// query_id = id->GetSeqIdString() + " " + kTitle;
/// if (query_id.size() > 35) query_id = query_id.substr(0, 25) + ".. ";
/// NCBI c++/src/algo/blast/api/blast_results.cpp:277-293:
/// retval = m_Errors.GetQueryId() + ": " + warning + " ";
///
/// `index` is the 0-based position of the query in the input; `program` is the lower-case
/// program name that NCBI's diagnostics print (`[tblastn]`, `[blastn]`).
pub fn invalid_query_warning(program: &str, index: usize, query: &fasta::Record) -> Vec<u8> {
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
    warning.extend_from_slice(b": Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options \n");
    warning
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
