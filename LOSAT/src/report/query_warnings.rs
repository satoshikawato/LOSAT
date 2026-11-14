//! Warnings that NCBI writes to standard error for a search, whatever the output format.

use crate::blastinput::fasta_reader::InputRecord;
use crate::report::archive::ArchiveHandle;
use crate::report::diag_log::Posts;
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
    before: Before<'a>,
    /// Standard error (the CLI) or the caller's diagnostics.
    pub sink: &'a mut (dyn Write + Send),
    /// The texts of the subjects that the reports read one by one (`with_subject_reads`), and
    /// how many of them are written.
    reads: Option<ReadCursor<'a>>,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
// ```c++
//        		const CBioseq_Handle& subject_bh =
//                 scope.GetBioseqHandle(align.GetSeq_id(1));
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1587-1589
// ```c++
//     try{
//         const CBioseq_Handle& handle = m_ScopeRef->GetBioseqHandle(*id);
//         x_FillDeflineAndId(handle, *id, use_this_seqid, sdl, blast_rank);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1947-1955
// ```c++
//             CRef<CAlnVec> avRef = x_GetAlnVecForSeqalign(**iter);
//     ...
//                 try{
//                     const CBioseq_Handle& handle = avRef->GetBioseqHandle(1);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:374-375
// ```c++
//     try {
//         const CBioseq_Handle& kSubjBioseqHandle = scope->GetBioseqHandle(kSeqId);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1465-1466
// ```c++
//     list< CRef<CSeqdesc> > tax =
//         x_GetTaxonomy(oid, target_gi, target_seq_id);
// ```
/// The taxonomy texts of the subjects that the queries' reports read one by one, where
/// `PRE_FETCH_SEQS_LIMIT` turns a batch's prefetch off: NCBI's formatter loads a subject when a
/// tabular row reads it, at its description line or, for one the description table does not
/// show, at its alignment's heading, and at its BLAST XML `Hit`; the first load of a blob looks
/// its taxids up (`x_GetTaxonomy`), and the posts of the lookup flush the report written before
/// it. By query index, the subjects whose first read writes a text, in the order the report
/// reads them, each with that text (made by the program's model of the reads,
/// `algorithm::common::database::subject_read_texts`).
#[derive(Debug, Default)]
pub struct SubjectReads {
    queries: Vec<QueryReads>,
}

/// The texts of one query's reads (`SubjectReads`): `(subject, text)` in the order of the
/// reads, a subject named as the report's writer names it (its index in the hit lists).
pub type QueryReads = Vec<(u32, Vec<u8>)>;

impl SubjectReads {
    /// Sets the texts of query `query`'s reads.
    pub fn set(&mut self, query: usize, reads: QueryReads) {
        if reads.is_empty() {
            return;
        }
        if self.queries.len() <= query {
            self.queries.resize_with(query + 1, Vec::new);
        }
        self.queries[query] = reads;
    }

    /// The texts of every query's reads, by query index (`QueryWarnings::with_subject_reads`).
    pub fn queries(&self) -> &[QueryReads] {
        &self.queries
    }
}

/// `reads` (by query) with each subject named by its rank in `ordinals` (the reported subjects'
/// ordinal ids, increasing): the name of the report writers of the programs that make the texts
/// before their hit lists name the subjects by rank (BLASTN, TBLASTX). A subject that is not
/// reported keeps no text.
pub fn reads_by_rank(reads: &[QueryReads], ordinals: &[u32]) -> Vec<QueryReads> {
    reads
        .iter()
        .map(|query| {
            query
                .iter()
                .filter_map(|(ordinal, text)| {
                    let rank = ordinals.binary_search(ordinal).ok()?;
                    Some((rank as u32, text.clone()))
                })
                .collect()
        })
        .collect()
}

/// The texts of query `query`'s reads in `reads`, in order.
fn reads_of(reads: &[QueryReads], query: usize) -> &[(u32, Vec<u8>)] {
    reads.get(query).map_or(&[], Vec::as_slice)
}

/// Every text of the reads of the queries `queries` in `reads`, in order.
pub fn read_texts(reads: &[QueryReads], queries: std::ops::Range<usize>) -> Vec<u8> {
    queries
        .flat_map(|query| reads_of(reads, query))
        .flat_map(|(_, text)| text.iter().copied())
        .collect()
}

/// A report writer's call where NCBI reads a subject (`QueryWarnings::subject_read`): the
/// subject's index and the writer, flushed before a text is written or, for a writer that
/// formats into a buffer, marked at the buffer's length.
pub type SubjectReadHook<'r, W> = &'r mut dyn FnMut(u32, &mut W) -> io::Result<()>;

/// A [`SubjectReadHook`] that writes nothing (the tests of the writers).
#[cfg(test)]
pub fn no_subject_reads<W>(_subject: u32, _writer: &mut W) -> io::Result<()> {
    Ok(())
}

/// The writing of the texts of `SubjectReads` where the reports read the subjects: every text
/// of the queries below `query` is written, and the first `next` of query `query`'s.
pub struct ReadCursor<'a> {
    reads: &'a [QueryReads],
    query: usize,
    next: usize,
}

impl<'a> ReadCursor<'a> {
    /// The cursor of `reads` (by query index), or `None` where no read writes a text.
    pub fn new(reads: &'a [QueryReads]) -> Option<Self> {
        reads
            .iter()
            .any(|query| !query.is_empty())
            .then_some(ReadCursor {
                reads,
                query: 0,
                next: 0,
            })
    }

    /// Whether the report of query `query` reads a subject that writes a text.
    pub fn has_texts(&self, query: usize) -> bool {
        !reads_of(self.reads, query).is_empty()
    }

    /// The report of query `query` reads subject `subject` here (a tabular row that reads it,
    /// its description line, its alignment's heading, its BLAST XML `Hit`): the texts of the
    /// query's reads up to its first are written to `sink`, after the report written so far
    /// (`report`) is flushed, unless they are written already. A subject whose reads write no
    /// text, or that the query's report read before, writes nothing.
    pub fn subject_read(
        &mut self,
        query: usize,
        subject: u32,
        report: &mut dyn Write,
        sink: &mut dyn Write,
    ) -> io::Result<()> {
        if query < self.query {
            return Ok(());
        }
        let from = if query == self.query { self.next } else { 0 };
        let reads = reads_of(self.reads, query);
        let Some(found) = reads
            .get(from..)
            .and_then(|rest| rest.iter().position(|(read, _)| *read == subject))
        else {
            return Ok(());
        };
        let mut text = self.unwritten_before(query);
        for (_, read) in &reads[from..=from + found] {
            text.extend_from_slice(read);
        }
        self.next = from + found + 1;
        write_after(report, sink, &text)
    }

    /// Writes to `sink` the texts of the reads of the queries below `query` that their reports
    /// did not place (`subject_read`), after the report written so far (`report`) is flushed:
    /// they come after those reports.
    pub fn before_query(
        &mut self,
        query: usize,
        report: &mut dyn Write,
        sink: &mut dyn Write,
    ) -> io::Result<()> {
        let text = self.unwritten_before(query);
        write_after(report, sink, &text)
    }

    /// The texts not written yet of every query (a report that ends without placing them).
    fn unwritten(&mut self) -> Vec<u8> {
        let queries = self.reads.len().max(self.query + 1);
        self.unwritten_before(queries)
    }

    /// The texts not written yet of the queries below `query`, which count as written after
    /// this; the cursor moves to the start of `query`'s.
    fn unwritten_before(&mut self, query: usize) -> Vec<u8> {
        if query <= self.query {
            return Vec::new();
        }
        let mut text = Vec::new();
        for (_, read) in reads_of(self.reads, self.query).iter().skip(self.next) {
            text.extend_from_slice(read);
        }
        text.extend(read_texts(self.reads, self.query + 1..query));
        (self.query, self.next) = (query, 0);
        text
    }
}

/// Writes the part `report` of query `query`'s report, formatted into a buffer, with the texts of
/// the subjects it reads one by one at its `reads` (the buffer's length where the serial writer
/// reads each, and the subject; `QueryWarnings::subject_read`): `write_hit_posts` of those marks
/// alone (no texts of the writer's own).
pub fn write_marked(
    report: &[u8],
    reads: &[(usize, u32)],
    query: usize,
    writer: &mut dyn Write,
    warnings: Option<&mut QueryWarnings<'_>>,
) -> io::Result<()> {
    let posts = HitPosts {
        every: true,
        marks: reads.to_vec(),
        ends: vec![0; reads.len()],
        texts: Vec::new(),
    };
    write_hit_posts(report, &posts, query, writer, warnings)
}

/// The places in a query's report, formatted into a buffer, where NCBI loads a subject and posts
/// texts: a structured report's `Hit` (BLAST XML, XML2, JSON; `HitSink::position`) or a tabular
/// row's or pairwise line's read (`write_marked`). At each, the text of the subject's read in the
/// program's model (`ReadCursor::subject_read`) and the texts that the report's writer makes there
/// (the taxonomy database's messages of the XML2 and JSON documents' `Hit`s, `TaxdbPrinter`).
/// NCBI posts them on `cerr`, which flushes the report written before them (`QueryWarnings`);
/// `write` writes them so.
#[derive(Debug, Default)]
pub struct HitPosts {
    /// Whether every `Hit` is marked (the model's reads of the query write texts); else only one
    /// whose writer makes a text.
    every: bool,
    /// By marked `Hit`: its place in the query's output and its subject.
    marks: Vec<(usize, u32)>,
    /// By marked `Hit`: the end of its texts in `texts`.
    ends: Vec<usize>,
    texts: Vec<u8>,
}

impl HitPosts {
    /// The posts of a query's `Hit`s; `every` where the model's reads of the query write texts
    /// (`QueryWarnings::reads_with_texts`).
    pub fn new(every: bool) -> Self {
        Self {
            every,
            ..Self::default()
        }
    }

    /// The `Hit` of subject `subject` loads it at `at` (`HitSink::position`); `post` writes the
    /// texts that the writer makes there.
    pub fn hit(
        &mut self,
        at: usize,
        subject: u32,
        post: impl FnOnce(&mut Vec<u8>) -> io::Result<()>,
    ) -> io::Result<()> {
        let from = self.texts.len();
        post(&mut self.texts)?;
        if self.every || self.texts.len() > from {
            self.marks.push((at, subject));
            self.ends.push(self.texts.len());
        }
        Ok(())
    }

    /// Writes `report`, the query's output, to `writer`, and at each marked `Hit`, after the
    /// bytes before it: `read` of its subject (the model's text, to `sink`; `writer` is the
    /// report) and the texts the writer made there, to `sink` after `writer` is flushed.
    pub fn write(
        &self,
        report: &[u8],
        writer: &mut dyn Write,
        sink: &mut dyn Write,
        mut read: impl FnMut(u32, &mut dyn Write, &mut dyn Write) -> io::Result<()>,
    ) -> io::Result<()> {
        let (mut index, mut from) = (0, 0);
        write_at_marks(report, &self.marks, writer, |subject, writer| {
            read(subject, writer, &mut *sink)?;
            let end = self.ends[index];
            index += 1;
            let texts = &self.texts[from..end];
            from = end;
            write_after(writer, &mut *sink, texts)
        })
    }
}

/// Writes query `query`'s structured report `report` with the posts of its `Hit`s
/// (`HitPosts::write`): the model's reads through the warnings' cursor
/// (`QueryWarnings::subject_read`) and the writer's texts, both to the warnings' sink. Without
/// warnings, the report alone.
pub fn write_hit_posts(
    report: &[u8],
    hits: &HitPosts,
    query: usize,
    writer: &mut dyn Write,
    warnings: Option<&mut QueryWarnings<'_>>,
) -> io::Result<()> {
    let Some(warnings) = warnings else {
        return writer.write_all(report);
    };
    let QueryWarnings { reads, sink, .. } = warnings;
    hits.write(
        report,
        writer,
        &mut **sink,
        |subject, writer, sink| match reads.as_mut() {
            Some(cursor) => cursor.subject_read(query, subject, writer, sink),
            None => Ok(()),
        },
    )
}

/// Writes `report`, a buffer, to `writer` with `read` of each subject at its place `marks` (the
/// buffer's length where the serial writer reads it, and the subject), after the bytes before it.
pub fn write_at_marks(
    report: &[u8],
    marks: &[(usize, u32)],
    writer: &mut dyn Write,
    mut read: impl FnMut(u32, &mut dyn Write) -> io::Result<()>,
) -> io::Result<()> {
    let mut written = 0;
    for &(at, subject) in marks {
        writer.write_all(&report[written..at])?;
        written = at;
        read(subject, writer)?;
    }
    writer.write_all(&report[written..])
}

/// Writes `text` to `sink` after `report` is flushed (nothing for an empty text).
fn write_after(report: &mut dyn Write, sink: &mut dyn Write, text: &[u8]) -> io::Result<()> {
    if text.is_empty() {
        return Ok(());
    }
    report.flush()?;
    sink.write_all(text)?;
    sink.flush()
}

impl Drop for QueryWarnings<'_> {
    /// The texts of the reads that the report did not place (a writer that stopped early, or a
    /// read no writer marks) are written when the warnings end, so that standard error keeps
    /// every text the reports' reads write.
    fn drop(&mut self) {
        if let Some(cursor) = self.reads.as_mut() {
            let text = cursor.unwritten();
            if !text.is_empty() {
                let _ = self.sink.write_all(&text);
                let _ = self.sink.flush();
            }
        }
    }
}

/// What `QueryWarnings` writes before the report of each query, by query index.
enum Before<'a> {
    /// Lines that are each a warning (`Posts::of_lines`): what the programs other than BLASTN
    /// post there.
    Lines(&'a [Vec<u8>]),
    /// The posts of each query, written with the BLAST archive's handle (`Posts::write`).
    Posts {
        posts: &'a [Posts],
        archive: Option<&'a ArchiveHandle>,
    },
}

impl<'a> QueryWarnings<'a> {
    /// The warnings `before` (lines, each a warning) of each query, written to `sink`.
    pub fn new(before: &'a [Vec<u8>], sink: &'a mut (dyn Write + Send)) -> Self {
        Self {
            before: Before::Lines(before),
            sink,
            reads: None,
        }
    }

    /// These warnings with the texts of the subjects that the reports read one by one
    /// (`SubjectReads`), which the report's writer places where it reads each subject
    /// (`subject_read`). A text that the writer does not place is written after the query's
    /// report, before the next query's warnings (`before_query`, `after_last_query`).
    pub fn with_subject_reads(mut self, reads: &'a [QueryReads]) -> Self {
        self.reads = ReadCursor::new(reads);
        self
    }

    /// Whether the reports of query `query` read a subject that writes a text (a writer that
    /// formats into buffers marks its reads only then, `subject_read`).
    pub fn reads_with_texts(&self, query: usize) -> bool {
        self.reads
            .as_ref()
            .is_some_and(|cursor| cursor.has_texts(query))
    }

    /// The texts of the subjects' reads (`with_subject_reads`; `None` without a read that writes
    /// one) and the sink, apart: a writer that writes texts of its own to the sink between the
    /// reads places them with the cursor (`tax_report::write_report`).
    pub fn reads_and_sink(&mut self) -> (Option<&mut ReadCursor<'a>>, &mut (dyn Write + Send)) {
        (self.reads.as_mut(), &mut *self.sink)
    }

    /// The report of query `query` reads subject `subject` here (`ReadCursor::subject_read`),
    /// after the report written so far (`report`).
    pub fn subject_read(
        &mut self,
        query: usize,
        subject: u32,
        report: &mut dyn Write,
    ) -> io::Result<()> {
        match self.reads.as_mut() {
            Some(cursor) => cursor.subject_read(query, subject, report, &mut *self.sink),
            None => Ok(()),
        }
    }

    /// The posts `before` of each query, written to `sink`, whose whole posts the BLAST archive's
    /// handle `archive` saves as one message each (`Posts::write`).
    pub fn of_posts(
        before: &'a [Posts],
        archive: Option<&'a ArchiveHandle>,
        sink: &'a mut (dyn Write + Send),
    ) -> Self {
        Self {
            before: Before::Posts {
                posts: before,
                archive,
            },
            sink,
            reads: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:304-308 and 562
    // ```c++
    //                         pssm =
    //                          ComputePssmForNextIteration(*seq, aln, psi_opts, scope,
    //                                      results_1st_query.GetAncillaryData());
    //     ...
    //         formatter.PrintEpilog(opt);
    // ```
    /// Writes the warnings that come after the report of the last of `queries` queries and before
    /// the epilog (psiblast's warnings of the last query's PSSM reads, which it makes after the
    /// query's report): `before[queries]` where the caller sets it, as `before_query`. A report
    /// writer calls it after its last query, before an epilog it writes (review 2 of line DBID
    /// pssmread, H4).
    pub fn after_last_query(&mut self, queries: usize, report: &mut dyn Write) -> io::Result<()> {
        self.before_query(queries, report)
    }

    /// Writes the warnings of query `index`, after the report written so far (`report`)
    /// is flushed; a failed flush returns before them.
    pub fn before_query(&mut self, index: usize, report: &mut dyn Write) -> io::Result<()> {
        // The texts of the reads that the reports of the queries before it did not place
        // (`subject_read`) come after those reports.
        if let Some(cursor) = self.reads.as_mut() {
            cursor.before_query(index, report, &mut *self.sink)?;
        }
        match &self.before {
            Before::Lines(before) => match before.get(index) {
                Some(warnings) if !warnings.is_empty() => {
                    report.flush()?;
                    self.sink.write_all(warnings)?;
                    self.sink.flush()
                }
                _ => Ok(()),
            },
            Before::Posts { posts, archive } => match posts.get(index) {
                Some(posts) if !posts.is_empty() => {
                    report.flush()?;
                    posts.write(&mut *self.sink, *archive)?;
                    self.sink.flush()
                }
                _ => Ok(()),
            },
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
#[cfg(test)]
pub fn invalid_query_warning<R: InputRecord>(program: &str, index: usize, query: &R) -> Vec<u8> {
    query_warning(program, index, query, &[INVALID_QUERY_MESSAGE.to_string()])
}

/// The warning of each query of the first query batch of a search with the IDENTITY matrix and a
/// composition mode above 0, which NCBI searches in mode 0 (blastp and tblastn; IDENTITY has no
/// composition data; the reset stays in the options of the later batches).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:150-154
/// ```c++
///         ostringstream os;
///         os << "Composition-based statistics cannot be used with the IDENTITY "
///            << "matrix, resetting the composition-based statistics option to 0";
///         search_messages.AddMessageAllQueries(eBlastSevWarning, 0, os.str());
/// ```
pub const IDENTITY_CBS_RESET_MESSAGE: &str = "Composition-based statistics cannot be used with the IDENTITY matrix, resetting the composition-based statistics option to 0";

/// NCBI's `kBlastErrMsg_CantCalculateUngappedKAParams` (see `invalid_query_warning`).
pub const INVALID_QUERY_MESSAGE: &str = "Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options";

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:207-214
// ```c++
//     if (! BlastSeqSrcGetNumSeqs(m_InternalData->m_SeqSrc->GetPointer())) {
//         string msg =
//             "Filtering resulted in an empty database.";
//
//         m_Messages.AddMessageAllQueries(eBlastSevWarning,
//                                         kBlastMessageNoContext,
//                                         msg);
//     }
// ```
/// The warning of each query of a searched batch when a list kept no sequence of the database
/// (`SeqDb::filtered_empty`). NCBI adds it in the preliminary search, after the messages of the
/// set-up are sorted, so it comes last (`sorted_query_messages`; probes C31 and D01 of line DI-8).
/// A batch without a valid context is not searched and does not have it.
pub const EMPTY_DATABASE_MESSAGE: &str = "Filtering resulted in an empty database.";

/// The messages of a query's result, which its warning writes to standard error
/// (`query_warning`) and the XML reports keep with the query (BLASTP and TBLASTN): the O
/// characters read as X (`replaced_o_message`), IDENTITY's reset of the composition mode
/// (`IDENTITY_CBS_RESET_MESSAGE`, for the queries of the first query batch), a record without
/// letters or else a query without a valid context (`INVALID_QUERY_MESSAGE`), and a database that
/// a list left empty (`EMPTY_DATABASE_MESSAGE`). `query_warning` sorts them as NCBI does.
pub fn query_messages(
    seq: &[u8],
    valid: bool,
    identity_cbs_reset: bool,
    empty_database: bool,
) -> Vec<String> {
    let mut messages = Vec::new();
    messages.extend(replaced_o_message(seq));
    if identity_cbs_reset {
        messages.push(IDENTITY_CBS_RESET_MESSAGE.to_string());
    }
    if seq.is_empty() {
        messages.push(crate::blastinput::query_batch::NO_DATA_MESSAGE.to_string());
    } else if !valid {
        messages.push(INVALID_QUERY_MESSAGE.to_string());
    }
    if empty_database {
        messages.push(EMPTY_DATABASE_MESSAGE.to_string());
    }
    messages
}

/// The messages of a query in the order of its warning: sorted without repeats (see
/// `query_warning`), `EMPTY_DATABASE_MESSAGE` last.
pub fn sorted_query_messages(messages: &[String]) -> Vec<&String> {
    let mut sorted: Vec<&String> = messages.iter().collect();
    sorted.sort_by_key(|message| (message.as_str() == EMPTY_DATABASE_MESSAGE, *message));
    sorted.dedup();
    sorted
}

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
    let local_id = format!("Query_{}", index + 1);
    let query_id = warning_query_id(
        query.believed_warning_id().unwrap_or(local_id.as_bytes()),
        &query.title_bytes(),
    );
    let mut warning = format!("Warning: [{program}] ").into_bytes();
    warning.extend_from_slice(&query_id);
    warning.extend_from_slice(b": ");
    for message in sorted_query_messages(messages) {
        warning.extend_from_slice(message.as_bytes());
        warning.push(b' ');
    }
    warning.push(b'\n');
    warning
}

/// `query_warning` for a query whose Seq-id is its own rather than the `Query_N` that NCBI's
/// FASTA reader makes: psiblast's PSSM query (`-in_pssm`), whose warnings start with the PSSM's
/// Seq-id, `seq_id` (`GetSeqIdString()`: a `local str` as written, a `local id` in decimal; the
/// PSSM path keeps it as the record's local id, `algorithm::psiblast::in_pssm`), then the title
/// `title` (the query source's `GetTitle`, which for a PSSM's Bioseq is not the record's title:
/// `blastp::blast_engine::pssm_hooks::pssm_query_warning`). The lines are `query_warning`'s, which
/// the other queries keep.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:533-535
/// ```c++
///             if (const CSeq_id* id = queries.GetSeqId(index)) {
///                 const string kTitle = queries.GetTitle(index);
///                 string query_id = id->GetSeqIdString();
/// ```
pub fn query_warning_with_seq_id(
    program: &str,
    seq_id: &[u8],
    title: &[u8],
    messages: &[String],
) -> Vec<u8> {
    if messages.is_empty() {
        return Vec::new();
    }
    let mut warning = format!("Warning: [{program}] ").into_bytes();
    warning.extend_from_slice(&warning_query_id(seq_id, title));
    warning.extend_from_slice(b": ");
    for message in sorted_query_messages(messages) {
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
    let list = residue_positions(&positions)?;
    Some(format!(
        "One or more O characters replaced by X for alignment score calculations at positions {list}"
    ))
}

/// The error of a protein query with residues that are no NCBIstdaa letter (`positions`: their
/// 0-based places in the converted sequence), or `None` without any. NCBI's `GetSequenceProtein`
/// throws it; `SetupQueries_OMF` catches it as a warning of the query and invalidates the query,
/// and when no query is valid the setup throws the messages of all queries
/// (`TSearchMessages::ToString`: each one with its severity, a space after it), which the
/// application reports as `BLAST engine error: <that text>` (exit 3).
///
/// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:906-918
/// ```c++
///     if (invalid_residues.size() > 0) {
///         string error("Invalid residues found at positions ");
///         error += NStr::IntToString(invalid_residues[0]);
///         for (i = 1; i < min(kMaxResiduesToWarnAbout, invalid_residues.size());
///              i++) {
///             error += ", " + NStr::IntToString(invalid_residues[i]);
///         }
///         if (invalid_residues.size() > kMaxResiduesToWarnAbout) {
///             error += ",... (only first ";
///             error += NStr::SizetToString(kMaxResiduesToWarnAbout) + " shown)";
///         }
///         NCBI_THROW(CBlastException, eInvalidCharacter, error);
///     }
/// ```
/// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:632-640 (the catch) and :651
/// (`NCBI_THROW(CBlastException, eSetup, messages.ToString())`), c++/include/algo/blast/api/blast_types.hpp:180-186
/// (`GetMessage(true)` is `GetSeverityString() + ": " + m_Message`).
pub fn invalid_residues_error(positions: &[usize]) -> Option<String> {
    let list = residue_positions(positions)?;
    Some(format!(
        "Warning: Invalid residues found at positions {list} "
    ))
}

/// The list of residue positions of NCBI's `GetSequenceProtein` messages: the first, then up to
/// 19 more after `, `, then `,... (only first 20 shown)` when there are more than 20
/// (`kMaxResiduesToWarnAbout`, blast_setup_cxx.cpp:920-932, quoted at `replaced_o_message`).
fn residue_positions(positions: &[usize]) -> Option<String> {
    let first = *positions.first()?;
    let mut list = first.to_string();
    for position in positions.iter().take(20).skip(1) {
        list.push_str(&format!(", {position}"));
    }
    if positions.len() > 20 {
        list.push_str(",... (only first 20 shown)");
    }
    Some(list)
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
    use crate::blastinput::fasta_reader::FastaRecord;

    #[test]
    fn messages_of_a_query_are_in_ncbis_sorted_order() {
        // NCBI blast_aux.cpp:1043-1054 sorts a query's messages (same error id and
        // severity, so by text) and removes repeats.
        let query = FastaRecord::new("Query_1", b"q1", b"OOO");
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

    use super::*;

    // NCBI 2.17.0 (M4a-1b probes s19, s20): the list form of `GetSequenceProtein`'s errors.
    #[test]
    fn the_invalid_residues_error_lists_the_first_twenty_positions() {
        assert_eq!(invalid_residues_error(&[]), None);
        assert_eq!(
            invalid_residues_error(&[26]).unwrap(),
            "Warning: Invalid residues found at positions 26 "
        );
        let twenty: Vec<usize> = (0..20).collect();
        let listed = (0..20)
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(", ");
        assert_eq!(
            invalid_residues_error(&twenty).unwrap(),
            format!("Warning: Invalid residues found at positions {listed} ")
        );
        let more: Vec<usize> = (0..22).collect();
        assert_eq!(
            invalid_residues_error(&more).unwrap(),
            format!(
                "Warning: Invalid residues found at positions {listed},... (only first 20 shown) "
            )
        );
    }

    // Probes C31 and D01 of line DI-8: NCBI adds the empty-database message to the sorted
    // messages of the set-up, though "Sequence ..." sorts after "Filtering ...".
    #[test]
    fn the_empty_database_message_comes_after_the_set_up_messages() {
        let messages = [
            EMPTY_DATABASE_MESSAGE.to_string(),
            "Sequence contains no data".to_string(),
            EMPTY_DATABASE_MESSAGE.to_string(),
        ];
        let query = FastaRecord::new("Query_2", b"empty", b"");
        assert_eq!(
            query_warning("blastp", 1, &query, &messages),
            b"Warning: [blastp] Query_2 empty: Sequence contains no data Filtering resulted in an empty database. \n"
        );
    }

    // NCBI 2.17.0 psiblast `-in_pssm` (the review of PS4a, S2 and N3; `PS/tmp/fix4a` probes
    // s2.idO, n3.localid_O, n3.localneg_O, n3.localid_notitle_O): the label is the PSSM's Seq-id,
    // whatever the query's place.
    #[test]
    fn a_pssm_querys_warning_starts_with_its_seq_id() {
        let o = replaced_o_message(b"MKVOAAOK").unwrap();
        let title = b"YP_009220601.1 TK-TMK, chimeric thymidine and thymidylate kinase";
        for (seq_id, title, shown) in [
            (
                &b"myq"[..],
                &title[..],
                &b"myq YP_009220601.1 TK-TMK.. "[..],
            ),
            (b"5", title, b"5 YP_009220601.1 TK-TMK, .. "),
            (b"-3", title, b"-3 YP_009220601.1 TK-TMK,.. "),
            (b"5", b"", b"5"),
        ] {
            let mut expected = b"Warning: [psiblast] ".to_vec();
            expected.extend_from_slice(shown);
            expected.extend_from_slice(b": ");
            expected.extend_from_slice(o.as_bytes());
            expected.extend_from_slice(b" \n");
            assert_eq!(
                query_warning_with_seq_id("psiblast", seq_id, title, std::slice::from_ref(&o)),
                expected
            );
        }
        assert!(query_warning_with_seq_id("psiblast", b"5", b"", &[]).is_empty());
    }

    // Query_1 + FASTA title is shortened after 35 bytes, then the invalid-Karlin
    // warning keeps NCBI's trailing space and newline.
    #[test]
    fn invalid_query_warning_matches_ncbi_bytes() {
        let short = FastaRecord::new("Query_1", b"nohit_query", b"W");
        assert_eq!(
            invalid_query_warning("tblastn", 0, &short),
            b"Warning: [tblastn] Query_1 nohit_query: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options \n"
        );
        let long = FastaRecord::new("Query_2", b"long_header abcdefghijklmnopqrstuvwxyz", b"W");
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

    /// The report (`R:`) and the warnings' texts in the order written.
    #[derive(Clone, Default)]
    struct Log(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);

    impl Write for Log {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(bytes);
            Ok(bytes.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn the_texts_of_the_reads_follow_the_report_up_to_each_read() {
        let reads: Vec<QueryReads> = vec![
            vec![
                (5, b"A\n".to_vec()),
                (7, b"B\n".to_vec()),
                (8, b"C\n".to_vec()),
            ],
            vec![(3, b"D\n".to_vec())],
            vec![(4, b"E\n".to_vec())],
        ];
        let before = vec![b"W0\n".to_vec(), b"W1\n".to_vec(), b"W2\n".to_vec()];
        let log = Log::default();
        let (mut report, mut sink) = (log.clone(), log.clone());
        let mut warnings = QueryWarnings::new(&before, &mut sink).with_subject_reads(&reads);
        assert!(warnings.reads_with_texts(0) && !warnings.reads_with_texts(3));
        warnings.before_query(0, &mut report).unwrap();
        // A subject without a text writes nothing; a read writes the texts up to its own.
        warnings.subject_read(0, 9, &mut report).unwrap();
        report.write_all(b"R1\n").unwrap();
        warnings.subject_read(0, 7, &mut report).unwrap();
        report.write_all(b"R2\n").unwrap();
        // Read before: nothing again.
        warnings.subject_read(0, 5, &mut report).unwrap();
        warnings.subject_read(0, 7, &mut report).unwrap();
        // C is never read: it comes after query 0's report, before query 1's warnings.
        warnings.before_query(1, &mut report).unwrap();
        report.write_all(b"R3\n").unwrap();
        warnings.subject_read(1, 3, &mut report).unwrap();
        // Query 2's text, never read, comes after the last report.
        warnings.before_query(2, &mut report).unwrap();
        warnings.after_last_query(3, &mut report).unwrap();
        drop(warnings);
        let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
        assert_eq!(text, "W0\nR1\nA\nB\nR2\nC\nW1\nR3\nD\nW2\nE\n");
    }

    #[test]
    fn the_posts_of_a_hit_follow_the_report_up_to_the_place_of_its_load() {
        // XML2 and JSON: every `Hit` loads its subject at the report's start, after outfmt 15's
        // separator; the report written before is flushed first (a buffered output stream).
        let log = Log::default();
        let mut report = io::BufWriter::with_capacity(1 << 16, log.clone());
        let mut sink = log.clone();
        report.write_all(b"{0}\n").unwrap();
        let mut posts = HitPosts::new(false);
        posts.hit(2, 1, |texts| texts.write_all(b"T1\n")).unwrap();
        posts.hit(2, 2, |_| Ok(())).unwrap();
        posts.hit(2, 3, |texts| texts.write_all(b"T3\n")).unwrap();
        posts
            .write(b",\n{1}\n", &mut report, &mut sink, |_, _, _| Ok(()))
            .unwrap();
        report.flush().unwrap();
        let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
        assert_eq!(text, "{0}\n,\nT1\nT3\n{1}\n");
        // BLAST XML: each `Hit` loads its subject at its own start; a query whose model's reads
        // write texts marks every `Hit`, and a `Hit`'s read comes before the texts made there.
        let reads: Vec<QueryReads> = vec![vec![(5, b"A\n".to_vec())]];
        let log = Log::default();
        let (mut report, mut sink) = (log.clone(), log.clone());
        let mut warnings = QueryWarnings::new(&[], &mut sink).with_subject_reads(&reads);
        let mut posts = HitPosts::new(warnings.reads_with_texts(0));
        posts.hit(0, 5, |texts| texts.write_all(b"P5\n")).unwrap();
        posts.hit(4, 6, |texts| texts.write_all(b"P6\n")).unwrap();
        write_hit_posts(b"<H5><H6>", &posts, 0, &mut report, Some(&mut warnings)).unwrap();
        drop(warnings);
        let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
        assert_eq!(text, "A\nP5\n<H5>P6\n<H6>");
        // Without warnings, the report alone.
        let mut alone = Vec::new();
        write_hit_posts(b"<H5><H6>", &posts, 0, &mut alone, None).unwrap();
        assert_eq!(alone, b"<H5><H6>");
    }

    #[test]
    fn a_query_without_texts_reads_nothing() {
        let reads: Vec<QueryReads> = vec![Vec::new(), Vec::new()];
        let mut sink = Vec::new();
        let warnings = QueryWarnings::new(&[], &mut sink).with_subject_reads(&reads);
        assert!(!warnings.reads_with_texts(0));
        assert!(ReadCursor::new(&reads).is_none());
    }
}
