//! Warnings that NCBI writes to standard error for a search, whatever the output format.

use crate::blastinput::fasta_reader::InputRecord;
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

/// Puts the reader's warnings about the query titles (`title_warning`) before the report of
/// the first query of each query batch. NCBI reads the queries one batch at a time, after
/// the outfmt 0 prolog and the reports of the previous batch, and its reader writes the
/// warnings of a batch's records as it reads them. `lengths` are the query lengths and
/// `batch_size` the residues of a batch (`query_batches`).
///
/// NCBI reference: c++/src/app/blast/blastp_app.cpp:253-259
/// ```c++
///         formatter.PrintProlog();
///
/// 	BLAST_PROF_ADD( BATCH_SIZE, (int)input.GetBatchSize() );
///         /*** Process the input ***/
///         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
/// 	    BLAST_PROF_START( APP.LOOP.PRE );
///             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
/// ```
pub fn prepend_batch_title_warnings<R: InputRecord, W: AsRef<[u8]>>(
    lines: &mut [Vec<u8>],
    queries: &[R],
    lengths: &[usize],
    batch_size: usize,
    title_warning: impl Fn(&R) -> W,
) {
    for range in crate::blastinput::query_batch::query_batches(lengths, batch_size) {
        let mut warnings: Vec<u8> = Vec::new();
        for query in &queries[range.clone()] {
            warnings.extend_from_slice(title_warning(query).as_ref());
        }
        if let Some(first) = lines.get_mut(range.start) {
            warnings.append(first);
            *first = warnings;
        }
    }
}

/// `prepend_batch_title_warnings` for the batches of a query input read with
/// `-query_loc` (`seq_range::QueryInput::batches`): NCBI's reader reads every record of a
/// batch, also those it skips, so their title warnings go before the first searched query
/// of the batch. `lines` are by searched query. Returns the title warnings of the last
/// batch when it has no searched query: NCBI writes them as it reads that batch, then
/// stops (`Empty CBlastQueryVector`).
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_input.cpp:144-155
/// ```c++
///         CRef<CBlastSearchQuery> q;
///         try { q.Reset(m_Source->GetNextSequence(scope)); }
///         ...
///         catch (const exception&) {
///             continue; //SB-2307. ignore well formed, not found accession
///         }
/// ```
pub fn prepend_ranged_batch_title_warnings<R: InputRecord, W: AsRef<[u8]>>(
    lines: &mut [Vec<u8>],
    input: &[R],
    batches: &[crate::blastinput::seq_range::RangedBatch],
    title_warning: impl Fn(&R) -> W,
) -> Vec<u8> {
    let mut unsearched = Vec::new();
    for batch in batches {
        let mut warnings: Vec<u8> = Vec::new();
        for query in &input[batch.input.clone()] {
            warnings.extend_from_slice(title_warning(query).as_ref());
        }
        match lines.get_mut(batch.searched.start) {
            Some(first) if !batch.searched.is_empty() => {
                warnings.append(first);
                *first = warnings;
            }
            _ => unsearched = warnings,
        }
    }
    unsearched
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
pub fn invalid_query_warning<R: InputRecord>(program: &str, index: usize, query: &R) -> Vec<u8> {
    query_warning(program, index, query, &[INVALID_QUERY_MESSAGE.to_string()])
}

/// NCBI's `kBlastErrMsg_CantCalculateUngappedKAParams` (see `invalid_query_warning`).
pub const INVALID_QUERY_MESSAGE: &str = "Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options";

/// The warning line of query `index` with its warning messages, each followed by a space
/// (`GetWarningStrings`, see `invalid_query_warning`); empty when there is none.
///
/// NCBI keeps a query's messages sorted and without repeats. Every message of one query
/// has the same severity (warning) and the same error id (the query's index), so NCBI's
/// `operator<` reduces to the comparison of the message strings.
///
/// NCBI reference: c++/src/algo/blast/api/blast_aux.cpp:1043-1054
/// ```c++
/// TSearchMessages::RemoveDuplicates()
/// {
///     NON_CONST_ITERATE(TSearchMessages, sm, (*this)) {
///         if (sm->empty()) {
///             continue;
///         }
///         sort(sm->begin(), sm->end(), TQueryMessagesLessComparator());
///         TQueryMessages::iterator new_end =
///             unique(sm->begin(), sm->end(), TQueryMessagesEqualComparator());
///         sm->erase(new_end, sm->end());
/// ```
/// NCBI reference: c++/include/algo/blast/api/blast_types.hpp:295-304
/// ```c++
/// CSearchMessage::operator<(const CSearchMessage& rhs) const
/// {
///     if (m_ErrorId < rhs.m_ErrorId ||
///         m_Severity < rhs.m_Severity ||
///         m_Message < rhs.m_Message) {
///         return true;
/// ```
/// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:620-622
/// ```c++
///                     CRef<CSearchMessage> m
///                         (new CSearchMessage(eBlastSevWarning, index, warnings));
///                     messages[index].push_back(m);
/// ```
pub fn query_warning<R: InputRecord>(
    program: &str,
    index: usize,
    query: &R,
    messages: &[String],
) -> Vec<u8> {
    if messages.is_empty() {
        return Vec::new();
    }
    let mut messages = messages.to_vec();
    messages.sort();
    messages.dedup();
    let query_id = warning_query_id(
        format!("Query_{}", index + 1).as_bytes(),
        &query.title_bytes(),
    );
    let mut warning = format!("Warning: [{program}] ").into_bytes();
    warning.extend_from_slice(&query_id);
    warning.extend_from_slice(b": ");
    for message in &messages {
        warning.extend_from_slice(message.as_bytes());
        warning.push(b' ');
    }
    warning.push(b'\n');
    warning
}

/// The query ID that NCBI's warnings of a query show: the local ID (`Query_N`), then a
/// space and the title when the query has one, cut to its first 25 bytes and `.. ` when
/// it is longer than 35 bytes. The bytes are the title's, raw; the cut can fall inside a
/// UTF-8 sequence.
///
/// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:533-543
/// ```c++
///             if (const CSeq_id* id = queries.GetSeqId(index)) {
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
pub fn warning_query_id(local_id: &[u8], title: &[u8]) -> Vec<u8> {
    let mut query_id = local_id.to_vec();
    if !title.is_empty() {
        query_id.push(b' ');
        query_id.extend_from_slice(title);
    }
    if query_id.len() > 35 {
        query_id.truncate(25);
        query_id.extend_from_slice(b".. ");
    }
    query_id
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
/// `query_warning` puts the query's messages in NCBI's order (its Karlin-Altschul message
/// sorts before this one).
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

    #[test]
    fn messages_of_a_query_are_in_ncbis_sorted_order() {
        // NCBI blast_aux.cpp:1043-1054 sorts a query's messages (same error id and
        // severity, so by text) and removes repeats.
        let query = bio::io::fasta::Record::with_attrs("q1", None, b"OOO");
        let o = replaced_o_message(b"OOO").unwrap();
        let line = query_warning(
            "blastp",
            0,
            &query,
            &[o.clone(), INVALID_QUERY_MESSAGE.to_string(), o.clone()],
        );
        let text = String::from_utf8(line).unwrap();
        assert!(
            text.starts_with("Warning: [blastp] Query_1 q1: Could not calculate"),
            "{text}"
        );
        assert_eq!(text.matches("One or more O").count(), 1, "{text}");
    }

    #[test]
    fn title_warnings_go_before_the_first_report_of_each_batch() {
        let records: Vec<_> = ["a", "b", "c"]
            .iter()
            .map(|id| bio::io::fasta::Record::with_attrs(id, None, b"ACDEFGHIKL"))
            .collect();
        let mut lines = vec![b"w0\n".to_vec(), Vec::new(), b"w2\n".to_vec()];
        // Batches of 15 residues: [a, b] and [c].
        prepend_batch_title_warnings(&mut lines, &records, &[10, 10, 10], 15, |record| {
            if record.id() == "c" {
                &b""[..]
            } else {
                b"T\n"
            }
        });
        assert_eq!(
            lines,
            vec![b"T\nT\nw0\n".to_vec(), Vec::new(), b"w2\n".to_vec()]
        );
    }

    use super::*;

    // Query_1 + FASTA title is shortened after 35 bytes, then the invalid-Karlin
    // warning keeps NCBI's trailing space and newline.
    #[test]
    fn invalid_query_warning_matches_ncbi_bytes() {
        let short = bio::io::fasta::Record::with_attrs("nohit_query", None, b"W");
        assert_eq!(
            invalid_query_warning("tblastn", 0, &short),
            b"Warning: [tblastn] Query_1 nohit_query: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options \n"
        );
        let long = bio::io::fasta::Record::with_attrs(
            "long_header",
            Some("abcdefghijklmnopqrstuvwxyz"),
            b"W",
        );
        assert!(invalid_query_warning("blastn", 1, &long)
            .starts_with(b"Warning: [blastn] Query_2 long_header abcde.. : "));
    }

    // NCBI BLAST+ 2.17.0, the invalid-query warning of BLASTN for one query with each
    // defline (inventory range RP, `scratch_RP/warn/iq_*_nuc_q.fa` and
    // `warn/out/iq_*_blastn_6.err`): the query ID between `[blastn] ` and `: Could`.
    #[test]
    fn warning_query_ids_follow_the_ncbi_oracle() {
        use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
        let w = |count: usize| vec![b'w'; count];
        let long = [&b"id "[..], &w(100)].concat();
        let cjk = [&[b'a'; 23][..], b"\xe6\xbc\xa2", &[b'z'; 14]].concat();
        let z27 = vec![b'z'; 27];
        let z28 = vec![b'z'; 28];
        for (defline, shown) in [
            (&b"sid sdesc"[..], &b"Query_1 sid sdesc"[..]),
            (b"", b"Query_1"),
            (b"id \xc3\xa9 x", b"Query_1 id \xc3\xa9 x"),
            (b"id \xe9 x", b"Query_1 id \xe9 x"),
            (b"id\tx", b"Query_1 id"),
            (&long, b"Query_1 id wwwwwwwwwwwwww.. "),
            (&cjk, b"Query_1 aaaaaaaaaaaaaaaaa.. "),
            (b"id x\xff", b"Query_1 id x\xff"),
            (&z27, b"Query_1 zzzzzzzzzzzzzzzzzzzzzzzzzzz"),
            (&z28, b"Query_1 zzzzzzzzzzzzzzzzz.. "),
        ] {
            let bytes = [&b">"[..], defline, b"\nACGTACGT\n"].concat();
            let mut source =
                FastaInputSource::from_bytes(&bytes, ReaderConfig::query("BLASTN", false, false));
            let records = read_all(&mut source, &mut |_| Ok(())).unwrap();
            assert_eq!(
                warning_query_id(records[0].local_id.as_bytes(), &records[0].title),
                shown,
                "{defline:?}"
            );
        }
    }
}
