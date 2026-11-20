//! The diagnostics that NCBI's application saves as the messages of a BLAST archive
//! (outfmt 11).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1155-1172
//! ```c++
//! void CBlastAppDiagHandler::Post(const SDiagMessage & mess)
//! {
//! 	if(m_handler != NULL) {
//! 		m_handler->Post(mess);
//! 	}
//! 	if(m_save) {
//! 		CRef<CBlast4_error> d(new CBlast4_error);
//! 		string m;
//! 		mess.Write(m);
//! 		d->SetMessage(NStr::Sanitize(m));
//! 		d->SetCode((int)mess.m_Severity);
//! 		{
//! 			DEFINE_STATIC_MUTEX(mx);
//! 			CMutexGuard guard(mx);
//! 			m_messages.push_back(d);
//! 		}
//! 	}
//! }
//! ```
//! Every diagnostic that the application posts is printed (the handler it wraps) and, when
//! the application writes an archive, saved with its severity (`eDiag_Warning` 1,
//! `eDiag_Error` 2, `eDiag_Critical` 3) and its text with the white space cleaned. LOSAT prints the same lines
//! (its standard error equals NCBI's), so the lines that go to its diagnostics stream are
//! the messages: `DiagLog` collects them as they are written, and an archive takes the ones
//! saved since the previous archive (`CBlastFormat::WriteArchive`, `m_Bah.ResetMessages()`).

/// `eDiag_Warning`.
pub const WARNING: i32 = 1;
/// `eDiag_Error`.
pub const ERROR: i32 = 2;
/// `eDiag_Critical`.
pub const CRITICAL: i32 = 3;

/// What NCBI posts at one place of a run, in order: standard error shows their texts one after
/// the other, and the BLAST archive saves each post as one message of its severity
/// (`CBlastAppDiagHandler::Post`). The diagnostics stream saves each line written to it as a
/// warning (`DiagTee`), which is what the reader's and the search's one-line warnings are: they
/// are kept as lines (`lines`). A post of another severity, or whose text has several lines (an
/// exception's `what()`), is kept whole (`post`) and saved as one message (`write`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Posts(Vec<Post>);

/// One item of `Posts`.
#[derive(Clone, Debug, PartialEq, Eq)]
enum Post {
    /// Lines that are each a warning.
    Lines(Vec<u8>),
    /// One post of `severity` (`WARNING`, `ERROR`, `CRITICAL`) with its text.
    Whole { severity: i32, text: Vec<u8> },
    /// The place of a search strategy request (`Posts::request`): writing it marks it reached.
    Request(RequestMark),
    /// The loads of a batch's prefetch (`Posts::loads`): lines whose text is made where they are
    /// written.
    Loads(Vec<crate::report::query_warnings::SubjectRead>),
    /// Lines of the query's result warnings that the call's first format does not write and a
    /// format after it does (`Posts::later`): they show nothing on the diagnostics.
    Later(Vec<u8>),
}

/// Where NCBI writes the search strategy request of the batch whose first query (the run's) is
/// `first` (`StrategyExport::write_with`), among the posts of the batch's start.
#[derive(Clone, Debug)]
struct RequestMark {
    export: crate::report::archive::strategy::StrategyExport,
    first: usize,
}

impl PartialEq for RequestMark {
    fn eq(&self, other: &Self) -> bool {
        self.first == other.first && self.export.same(&other.export)
    }
}

impl Eq for RequestMark {}

impl Posts {
    /// No post.
    pub const fn new() -> Self {
        Self(Vec::new())
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-282
    // ```c++
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    //             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
    //
    //             SaveSearchStrategy(args, m_CmdLineArgs, queries, m_OptsHndl);
    // ```
    /// The place where NCBI writes the search strategy request of `export`'s batch whose first
    /// query (the run's) is `first`: writing these posts there (`write`) marks the requests through
    /// that batch as NCBI's (`StrategyExport::release_through`; those after the last place reached
    /// are cut from the file when NCBI's error ends the run). It shows nothing on the diagnostics
    /// (`has_text`).
    pub fn request(export: crate::report::archive::strategy::StrategyExport, first: usize) -> Self {
        Self(vec![Post::Request(RequestMark { export, first })])
    }

    /// The posts of `lines`, each line a warning.
    pub fn of_lines(lines: Vec<u8>) -> Self {
        let mut posts = Self::new();
        if !lines.is_empty() {
            posts.0.push(Post::Lines(lines));
        }
        posts
    }

    /// Adds `lines`, each line a warning.
    pub fn lines(&mut self, lines: &[u8]) {
        if lines.is_empty() {
            return;
        }
        match self.0.last_mut() {
            Some(Post::Lines(last)) => last.extend_from_slice(lines),
            _ => self.0.push(Post::Lines(lines.to_vec())),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1432-1452
    // ```c++
    //     if (m_FormatType == CFormattingArgs::eAsnText
    //       || m_FormatType == CFormattingArgs::eAsnBinary
    //       || m_FormatType == CFormattingArgs::eXml
    //     ...
    //       || m_FormatType == CFormattingArgs::eSAM)
    //     {
    //         x_PrintStructuredReport(results, queries);
    //         return;
    //     }
    //     ...
    //     if (results.HasWarnings()) {
    //         ERR_POST(Warning << results.GetWarningStrings());
    //     }
    // ```
    /// Adds `lines`, the query's result warnings, where the call's first format is a structured
    /// one, which returns before them, and the caller takes those of the formats after it
    /// (`ReportOutputs::later_warnings`): the command line's run of such a format posts them when
    /// the query's report starts, so its writer writes them to the format's receiver there
    /// (`write_later`, `QueryWarnings::with_later`). Nothing of them shows on the diagnostics
    /// (`write`, `text`, `has_text`), which are the first format's.
    pub fn later(&mut self, lines: &[u8]) {
        if !lines.is_empty() {
            self.0.push(Post::Later(lines.to_vec()));
        }
    }

    /// `posts`, the posts before each query's report, with the later lines of each query
    /// (`later`, with NCBI's lines above it; `lines` by query, none for a query beyond it).
    pub fn with_later(mut posts: Vec<Posts>, lines: &[Vec<u8>]) -> Vec<Posts> {
        for (posts, lines) in posts.iter_mut().zip(lines) {
            posts.later(lines);
        }
        posts
    }

    /// Writes the lines of `later` (with NCBI's lines above it) to `receiver`, a later format's
    /// (`ReportOutputs::later_warnings`).
    pub fn write_later(&self, receiver: &mut Vec<u8>) {
        for post in &self.0 {
            if let Post::Later(lines) = post {
                receiver.extend_from_slice(lines);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288
    // ```c++
    //                 BlastFormatter_PreFetchSequenceData(*results, scope,
    //                 		                            fmt_args->GetFormattedOutputChoice());
    //                 ITERATE(CSearchResultSet, result, *results) {
    //                     formatter.PrintOneResultSet(**result, query_batch);
    //                 }
    // ```
    /// Adds the texts of the loads of a batch's prefetch (`reads`, in the order of the loads;
    /// `common::database::batch_prefetch`), each line a warning. NCBI's prefetch loads the batch's
    /// subjects after the reports of the batches before, whose loads a BLAST XML stream that cut
    /// a report did not make (`query_warnings::ReadCursor::cut`): their texts go to the prefetch's
    /// loads that make the same asks, so each load's text is made where the posts are written
    /// (`SubjectRead::placed`), not where the program models the batch.
    pub fn loads(&mut self, reads: Vec<crate::report::query_warnings::SubjectRead>) {
        if reads.iter().any(|read| read.load.is_some()) {
            self.0.push(Post::Loads(reads));
        } else {
            // No text is left to a load of a run whose lookups write none (`ReadLoad`).
            let text: Vec<u8> = reads.into_iter().flat_map(|read| read.text).collect();
            self.lines(&text);
        }
    }

    /// The reads of the loads of a batch's prefetch among these posts (`loads`), in order.
    pub fn load_reads(&self) -> impl Iterator<Item = &crate::report::query_warnings::SubjectRead> {
        self.0.iter().flat_map(|post| match post {
            Post::Loads(reads) => reads.as_slice(),
            _ => &[],
        })
    }

    /// Adds one post of `severity` whose text is `text` (with its line feed).
    pub fn post(&mut self, severity: i32, text: impl Into<Vec<u8>>) {
        let text = text.into();
        if !text.is_empty() {
            self.0.push(Post::Whole { severity, text });
        }
    }

    /// Adds the posts of `other` after these.
    pub fn append(&mut self, other: &Posts) {
        for post in &other.0 {
            match post {
                Post::Lines(lines) => self.lines(lines),
                Post::Whole { .. } | Post::Request(_) | Post::Loads(_) | Post::Later(_) => {
                    self.0.push(post.clone())
                }
            }
        }
    }

    /// These posts with the texts of their loads (`loads`) made here, once: the place where they
    /// are written. The texts that the loads take from the run's loads left out are theirs from
    /// now on.
    pub fn placed(&self) -> Posts {
        self.placed_in(None)
    }

    /// [`Posts::placed`] with the texts left out in `ledger` where it is given instead of the
    /// run's (`query_warnings::SubjectRead::placed_in`: the model of the writers).
    pub fn placed_in(
        &self,
        ledger: Option<&std::sync::Arc<crate::report::query_warnings::LoadTexts>>,
    ) -> Posts {
        let mut posts = Posts::new();
        for post in &self.0 {
            match post {
                Post::Loads(reads) => {
                    let text: Vec<u8> = reads
                        .iter()
                        .flat_map(|read| read.placed_in(ledger))
                        .collect();
                    posts.lines(&text);
                }
                Post::Lines(lines) => posts.lines(lines),
                Post::Whole { .. } | Post::Request(_) | Post::Later(_) => {
                    posts.0.push(post.clone())
                }
            }
        }
        posts
    }

    /// Whether there is no post (and no request mark, `request`).
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Whether a post shows text on the diagnostics (a request mark alone shows none, nor do a
    /// later format's lines, `later`; loads, whose texts are made where they are written, count
    /// as text unless none of them can make one: `placed` first where it decides).
    pub fn has_text(&self) -> bool {
        self.0.iter().any(|post| match post {
            Post::Request(_) | Post::Later(_) => false,
            Post::Loads(reads) => reads
                .iter()
                .any(|read| !read.text.is_empty() || read.load.is_some()),
            Post::Lines(_) | Post::Whole { .. } => true,
        })
    }

    /// What standard error shows of the posts.
    pub fn text(&self) -> Vec<u8> {
        let mut text = Vec::new();
        for post in &self.0 {
            match post {
                Post::Lines(bytes) | Post::Whole { text: bytes, .. } => {
                    text.extend_from_slice(bytes)
                }
                Post::Loads(reads) => text.extend(reads.iter().flat_map(|read| read.placed())),
                Post::Request(_) | Post::Later(_) => {}
            }
        }
        text
    }

    /// Writes the posts to `diagnostics`, the stream of the run's diagnostics: the lines as they
    /// are, and with the BLAST archive (`archive`, whose handle saves what that stream writes)
    /// each whole post saved as one message (`ArchiveHandle::write_post`); a request mark marks its
    /// search strategy requests reached (`request`); the loads' texts are made here (`placed`). A
    /// later format's lines (`later`) are not the diagnostics': nothing of them is written.
    pub fn write(
        &self,
        diagnostics: &mut dyn std::io::Write,
        archive: Option<&crate::report::archive::ArchiveHandle>,
    ) -> std::io::Result<()> {
        for post in &self.0 {
            match (post, archive) {
                (Post::Request(mark), _) => mark.export.release_through(mark.first)?,
                (Post::Later(_), _) => {}
                (Post::Loads(reads), _) => {
                    for read in reads {
                        diagnostics.write_all(&read.placed())?;
                    }
                }
                (Post::Whole { severity, text }, Some(handle)) => {
                    handle.write_post(diagnostics, *severity, text)?
                }
                (Post::Lines(text) | Post::Whole { text, .. }, _) => diagnostics.write_all(text)?,
            }
        }
        Ok(())
    }
}

/// A saved message: its severity and its sanitized text.
pub type Message = (i32, Vec<u8>);

/// The messages saved since the last archive, and the line being written.
#[derive(Clone, Debug, Default)]
pub struct DiagLog {
    saved: Vec<Message>,
    line: Vec<u8>,
}

impl DiagLog {
    /// Saves the lines that `bytes` complete as warnings (a line of standard error is a
    /// warning or reader message; the error lines of an aborted run are posted with
    /// `post`).
    pub fn feed(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            if byte == b'\n' {
                let line = std::mem::take(&mut self.line);
                self.post(WARNING, &line);
            } else {
                self.line.push(byte);
            }
        }
    }

    /// Saves a message of the given severity (`NStr::Sanitize` of its text).
    pub fn post(&mut self, severity: i32, line: &[u8]) {
        let text = sanitize(line);
        if !text.is_empty() {
            self.saved.push((severity, text));
        }
    }

    /// The messages saved since the last call (`ResetMessages` after the archive).
    pub fn take(&mut self) -> Vec<Message> {
        std::mem::take(&mut self.saved)
    }

    /// Whether no message is waiting for an archive.
    pub fn is_empty(&self) -> bool {
        self.saved.is_empty()
    }
}

/// `NStr::Sanitize(str)` with its default flags (`fSS_print`): a printable ASCII byte stays,
/// any other byte becomes a space, runs of spaces become one, and the spaces at both ends
/// go.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3605-3718
/// (`isprint` of the C locale is 0x20..=0x7e; a rejected byte is replaced by the space
/// `reject_replacement`, merged with the space before it, and leading and trailing spaces
/// are dropped).
pub fn sanitize(line: &[u8]) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::with_capacity(line.len());
    for &byte in line {
        let byte = if (0x20..=0x7e).contains(&byte) {
            byte
        } else {
            b' '
        };
        if byte == b' ' && (out.is_empty() || out.last() == Some(&b' ')) {
            continue;
        }
        out.push(byte);
    }
    if out.last() == Some(&b' ') {
        out.pop();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sanitize_cleans_white_space_like_ncbi() {
        assert_eq!(sanitize(b"  a  b\t\tc \r"), b"a b c");
        assert_eq!(sanitize(b"x\x01\x02y"), b"x y");
        assert_eq!(sanitize(b"caf\xc3\xa9!"), b"caf !");
        assert_eq!(sanitize(b" \t "), b"");
        assert_eq!(
            sanitize(b"File is not accessible:  `a b'"),
            b"File is not accessible: `a b'"
        );
        assert_eq!(sanitize(b"\"quoted\""), b"\"quoted\"");
    }

    #[test]
    fn lines_are_saved_when_complete() {
        let mut log = DiagLog::default();
        log.feed(b"Warning: [blastn] one\nWar");
        log.feed(b"ning: [blastn]  two\n\n");
        log.post(ERROR, b"BLAST query error: x\n");
        assert_eq!(
            log.take(),
            vec![
                (WARNING, b"Warning: [blastn] one".to_vec()),
                (WARNING, b"Warning: [blastn] two".to_vec()),
                (ERROR, b"BLAST query error: x".to_vec()),
            ]
        );
        assert!(log.is_empty());
    }

    // The later formats' lines of a query (`Posts::later`) are not the diagnostics': they are
    // written to a later format's receiver alone, before the query's report
    // (`QueryWarnings::with_later`), and kept by the copies of the posts.
    #[test]
    fn later_lines_go_to_a_later_format_s_receiver_alone() {
        use crate::report::query_warnings::QueryWarnings;
        let mut first = Posts::of_lines(b"Warning: [blastn] reader\n".to_vec());
        first.later(b"Warning: [blastn] Query_1 q: Sequence contains no data \n");
        let mut second = Posts::new();
        second.later(b"");
        assert!(second.is_empty());
        second.later(b"Warning: [blastn] Query_2 r: Sequence contains no data \n");
        assert!(!second.has_text() && second.text().is_empty());
        let posts = Posts::with_later(vec![first, second, Posts::new()], &[Vec::new()]);
        let copies: Vec<Posts> = posts.iter().map(Posts::placed).collect();
        assert_eq!(copies, posts);
        // The first format: the diagnostics get the reader's line, no receiver anything.
        let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
        let mut warnings = QueryWarnings::of_posts(&posts, None, &mut diagnostics);
        for query in 0..3 {
            warnings.before_query(query, &mut report).unwrap();
        }
        drop(warnings);
        assert_eq!(diagnostics, b"Warning: [blastn] reader\n");
        // A later format whose report stops after the first query: the diagnostics get nothing,
        // its receiver the first query's line.
        let (mut diagnostics, mut later) = (Vec::new(), Vec::new());
        let mut warnings = QueryWarnings::of_posts(&[], None, &mut diagnostics)
            .with_later(&posts, Some(&mut later));
        warnings.before_query(0, &mut report).unwrap();
        drop(warnings);
        assert!(diagnostics.is_empty() && report.is_empty());
        assert_eq!(
            later,
            b"Warning: [blastn] Query_1 q: Sequence contains no data \n"
        );
        // One that writes every query's report: both lines, in the order of the queries.
        let mut later = Vec::new();
        let mut warnings = QueryWarnings::of_posts(&[], None, &mut diagnostics)
            .with_later(&posts, Some(&mut later));
        for query in 0..3 {
            warnings.before_query(query, &mut report).unwrap();
        }
        drop(warnings);
        assert_eq!(
            later,
            b"Warning: [blastn] Query_1 q: Sequence contains no data \n\
              Warning: [blastn] Query_2 r: Sequence contains no data \n"
        );
    }
}
