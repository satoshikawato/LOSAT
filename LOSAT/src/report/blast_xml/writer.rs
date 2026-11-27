//! The type writers of the `CBlastOutput` model (objects/blastxml/blastxml.asn) over an
//! `Encoder`, and BLAST's incremental document (`BlastXML_FormatReport`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:725-733
//! So a document is: once, the `BlastOutput` object (filled from the first formatted query)
//! cut after `</BlastOutput_param>` and `\n<BlastOutput_iterations>`; per query the
//! `Iteration` object cut after `</Iteration_query-len>`, `\n<Iteration_hits>\n`, each `Hit`
//! written alone, `</Iteration_hits>` and the `Iteration` tail (`Statistics`, `message`,
//! `\n</Iteration>`); at the end `\n</BlastOutput_iterations>\n</BlastOutput>\n` and the
//! epilog's end of line (`CBlastFormat::PrintEpilog`, blast_format.cpp:2242-2243).

use anyhow::Result;

use crate::report::serial::encoder::Encoder;
use crate::report::serial::seqalign::AlignRow;
use crate::report::serial::xml::{XmlEncoder, XmlPost, XmlThrow};

use super::hsp::{mask_query_row, recount, round6, MaskLetter};
use super::names::NO_DEFINITION_LINE;
use super::records::{QueryMask, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats};
use crate::common::GapEditOp;

pub(crate) fn string_member<E: Encoder>(enc: &mut E, name: &str, tag: u8, value: &[u8]) {
    enc.member(name, tag);
    enc.visible_str(value);
    enc.end_member();
}

pub(crate) fn int_member<E: Encoder>(enc: &mut E, name: &str, tag: u8, value: i32) {
    enc.member(name, tag);
    enc.int(value);
    enc.end_member();
}

pub(crate) fn real_member<E: Encoder>(enc: &mut E, name: &str, tag: u8, value: f64) {
    enc.member(name, tag);
    enc.real(value);
    enc.end_member();
}

/// `Parameters` (`s_SetBlastXMLParameters`, blastxml_format.cpp:598-626).
fn write_parameters<E: Encoder>(enc: &mut E, run: &XmlRun<'_>) {
    enc.member("param", 8);
    enc.begin_sequence("Parameters");
    if let Some(matrix) = run.matrix {
        string_member(enc, "matrix", 0, matrix.as_bytes());
    }
    real_member(enc, "expect", 1, run.expect);
    if let Some((reward, penalty)) = run.reward_penalty {
        // NCBI writes sc-match only when `GetMatchReward` is not 0.
        if reward != 0 {
            int_member(enc, "sc-match", 3, reward);
        }
        if penalty != 0 {
            int_member(enc, "sc-mismatch", 4, penalty);
        }
    }
    int_member(enc, "gap-open", 5, run.gap_open);
    int_member(enc, "gap-extend", 6, run.gap_extend);
    if !run.filter.is_empty() {
        string_member(enc, "filter", 7, run.filter);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:616-618
    if let Some(pattern) = run.pattern {
        string_member(enc, "pattern", 8, pattern);
    }
    enc.end_sequence();
    enc.end_member();
}

/// What NCBI's XML writers post on the diagnostics stream while they write a document (the 803.7
/// post of a replaced control character, the report of an exception they catch): `text`, after the
/// bytes of the call's output buffer before `at` (`StructuredXml::take_posts`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DocumentPost {
    pub at: usize,
    pub text: Vec<u8>,
}

/// The 803.7 posts of an object's replaced control characters (`XmlPost`), each after `at(post)`
/// bytes of the output; `program` is the application's prefix of the posts.
pub(crate) fn replaced_posts(
    posts: Vec<XmlPost>,
    program: &str,
    mut at: impl FnMut(&XmlPost) -> usize,
) -> Vec<DocumentPost> {
    posts
        .into_iter()
        .map(|post| DocumentPost {
            at: at(&post),
            text: post.replaced.post(Some(program)),
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
/// An exception of NCBI's XML stream, or of the making of a 13 or 14 side file's report, that ends
/// the run (`CATCH_ALL`: `Error: `, `what()`, exit status 255; a `CSeqDBException` of a subject's
/// lookup `BLAST Database error: `, its message, exit status 2; `cli::NativeError::exception`).
/// What NCBI posted before it (the 803.7 posts of the object, the texts of the loads before it)
/// is written to the run's diagnostics first, by the one writer of a document's posts
/// (`report::query_warnings::write_hit_posts`), which hands this error back after them
/// (`DocumentPosts::ends_run`).
pub(crate) fn run_ending_exception(
    exception: &crate::ncbi_exception::NcbiException,
) -> anyhow::Error {
    crate::cli::NativeError::exception(exception).into()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:768-778,877-886 (cited at
// `blast_xml2::Xml2Document::failed`), 877-888
/// The outcome `written` of one query's structured report (`StructuredXml::query` of
/// `document`): `Ok(None)` where the report is written; where it raised an exception that NCBI's
/// report function catches (`CaughtException`: the stream forms of the XML2 and JSON documents,
/// outfmt 16 and 15, leave the report out, `blast_xml2::Xml2Document::report`), `Ok(Some)` with
/// the post of it (`ERR_POST` at the Error level after the prefix `[app] `, the running
/// application's name `app`,
/// `StructuredXml::caught_message`), which comes after the posts of the loads that the report made
/// before it (`query_warnings::HitPosts::caught`); any other error is returned, and the caller
/// drops what the failed report left in its buffer.
pub fn caught_post<D: StructuredXml + ?Sized>(
    written: Result<()>,
    document: &D,
    app: &str,
) -> Result<Option<Vec<u8>>> {
    match written {
        Ok(()) => Ok(None),
        Err(error) => match error.downcast::<super::CaughtException>() {
            Ok(super::CaughtException(exception)) => {
                Ok(Some(crate::ncbi_exception::error_post_line(
                    None,
                    Some(app),
                    &document.caught_message(&exception),
                )))
            }
            Err(error) => Err(error),
        },
    }
}

/// What a document's writer hands its caller with the bytes it wrote since the last
/// `StructuredXml::take_posts`: its posts (`DocumentPost`, at offsets of the output buffers of the
/// calls that wrote them), whether a query's report left out the loads of its later subjects, and
/// an exception that ends the run once the posts are written. The caller writes them through
/// `report::query_warnings::write_hit_posts` (or `write_document_posts`), the one place that
/// writes a document's posts and hands the run's end back. Dropping it unwritten loses the
/// posts and NCBI's end of the run, so it is `must_use`.
#[derive(Debug, Default)]
#[must_use = "a document's posts and the exception that ends the run go through query_warnings::write_hit_posts"]
pub struct DocumentPosts {
    pub posts: Vec<DocumentPost>,
    /// The query written last made no `Hit` after the last one marked (`HitPosts::hit`), so NCBI
    /// loaded none of their subjects (`s_SeqAlignToXMLHit` loads a subject when it makes its
    /// `Hit`): in BLAST XML (outfmt 5) its `Iteration` threw, or one of its `Hit`s did; in BLAST
    /// XML, XML2 and JSON one of its HSPs raised an exception that the report function catches
    /// (`CaughtException`).
    pub cut: bool,
    /// The exception that ends the run after these posts: of NCBI's XML stream (the BLAST XML
    /// header, a 14 side file) or of the making of a 13 or 14 side file's report
    /// (`Xml2Document::failed`) (`run_ending_exception`).
    pub ends_run: Option<Box<crate::ncbi_exception::NcbiException>>,
    /// Where the post of a held report's exception that the stream forms catch goes in the
    /// epilog's bytes (`Xml2Document::epilog`): after the report's separator, before the
    /// document's end, which NCBI writes after the post (`place_caught`).
    pub caught_at: Option<usize>,
    /// Where NCBI's insertions into the report's stream end, and its regions and flushes, at
    /// offsets of the output buffers (`ostream::Cuts`; BLAST XML's, `BlastXmlWriter`).
    pub cuts: crate::report::ostream::Cuts,
}

impl DocumentPosts {
    /// Puts `post`, the post of a held report's exception that the stream forms catch, where the
    /// epilog marked it (`caught_at`: before the document's end); `Some(post)` back where nothing
    /// marked it, for the caller to write after the report (`query_warnings::HitPosts::caught`).
    pub fn place_caught(&mut self, post: Vec<u8>) -> Option<Vec<u8>> {
        let Some(at) = self.caught_at.take() else {
            return Some(post);
        };
        self.posts.push(DocumentPost { at, text: post });
        self.posts.sort_by_key(|post| post.at);
        None
    }
}

/// The marker error with which outfmt 5's hit level stops the program's hits of a query whose
/// `Hit` threw: NCBI's exception leaves `s_SeqAlignSetToXMLHits` there, so no later hit is made
/// (`BlastXmlWriter::query` catches it).
#[derive(Debug)]
struct IterationThrew;

impl std::fmt::Display for IterationThrew {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("the XML iteration threw")
    }
}

impl std::error::Error for IterationThrew {}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:740-764
/// The post of an exception in an `Iteration` or one of its `Hit`s: each query is formatted by
/// itself (`CCmdLineBlastXMLReportData` of one result), so the label is always `Q(0/1` (sic).
fn iteration_post(at: usize, app: &str, thrown: &XmlThrow) -> DocumentPost {
    let message = thrown
        .exception()
        .what_after(b"Failed s_BlastXMLAddIteration Q(0/1");
    DocumentPost {
        at,
        text: crate::ncbi_exception::error_post_line(None, Some(app), &message),
    }
}

/// The `BlastOutput` object through `</BlastOutput_param>`, filled from the first formatted
/// query (`BlastXML_FormatReport`, blastxml_format.cpp:679-720), then `\n<BlastOutput_iterations>`.
/// NCBI serializes it into a string first, outside any `try` (c++/src/algo/blast/format/blastxml_format.cpp:724-733): its
/// posts come before it, and an exception ends the run with nothing of it written: the exception
/// is returned (`Ok(Some)`, `DocumentPosts::ends_run`), after the posts.
fn write_header(
    out: &mut Vec<u8>,
    run: &XmlRun<'_>,
    first: &XmlQuery<'_>,
    posts: &mut Vec<DocumentPost>,
) -> Result<Option<Box<XmlThrow>>> {
    let start = out.len();
    let mut enc = XmlEncoder::document(out, "NCBI-BlastOutput");
    enc.begin_document("BlastOutput");
    enc.begin_sequence("BlastOutput");
    string_member(&mut enc, "program", 0, run.program.as_bytes());
    string_member(&mut enc, "version", 1, run.version.as_bytes());
    string_member(&mut enc, "reference", 2, run.reference.as_bytes());
    string_member(&mut enc, "db", 3, run.db_name);
    string_member(&mut enc, "query-ID", 4, first.id);
    string_member(&mut enc, "query-def", 5, query_def(first.title));
    int_member(&mut enc, "query-len", 6, first.len);
    write_parameters(&mut enc, run);
    enc.raw(b"\n<BlastOutput_iterations>");
    let header_posts = enc.take_posts();
    let thrown = enc.take_thrown();
    let error = enc.error();
    posts.extend(replaced_posts(header_posts, run.app, |_| start));
    if let Some(error) = error {
        out.truncate(start);
        return match thrown {
            Some(thrown) => Ok(Some(thrown)),
            None => Err(error.into()),
        };
    }
    Ok(None)
}

/// `Statistics` (`s_BlastXMLGetStatistics`, blastxml_format.cpp:628-650).
fn write_statistics<E: Encoder>(enc: &mut E, run: &XmlRun<'_>, stats: &XmlStats) {
    enc.member("stat", 5);
    enc.begin_sequence("Statistics");
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:116-119
    // `m_NumSequences` is an `int`: the `Int8` sum fits it (`XmlRun::db_num`).
    int_member(enc, "db-num", 0, run.db_num as i32);
    enc.member("db-len", 1);
    enc.big_int(run.db_len);
    enc.end_member();
    int_member(enc, "hsp-len", 2, stats.hsp_len);
    // `stat->SetEff_space` takes `GetEffectiveSearchSpace` of `index`, cast to a double
    real_member(enc, "eff-space", 3, stats.eff_space as f64);
    real_member(enc, "kappa", 4, stats.kappa);
    real_member(enc, "lambda", 5, stats.lambda);
    real_member(enc, "entropy", 6, stats.entropy);
    enc.end_sequence();
    enc.end_member();
}

/// `Iteration_query-def` / `Hit_def`: the generated title, or `No definition line` for an
/// empty one (blastxml_format.cpp:557-565).
fn query_def(title: &[u8]) -> &[u8] {
    if title.is_empty() {
        NO_DEFINITION_LINE
    } else {
        title
    }
}

/// Buffers reused from HSP to HSP and query to query.
#[derive(Default)]
pub(crate) struct Scratch {
    pub(crate) masked: Vec<u8>,
    /// The subject row with its end gaps as spaces, where it has some (`HspCounts::hseq`).
    pub(crate) hseq: Vec<u8>,
    pub(crate) midline: Vec<u8>,
    pub(crate) columns: Vec<usize>,
    pub(crate) segments: Vec<(usize, usize)>,
    hit: Vec<u8>,
    /// The `Iteration`'s part after `</Iteration_query-len>` (`BlastXmlWriter::query`).
    tail: Vec<u8>,
}

/// The writer of a structured BLAST report (outfmt 5, and the XML2 / JSON formats 13-16): a
/// program's adapter hands it each query and its hits one at a time through the sinks, and
/// finishes the document. Each program's adapter is written once over it.
pub trait StructuredXml {
    /// Appends the report of one query to `out` (before it the document's start, for the
    /// first query). `hits` writes the query's hits, one subject at a time, through the sink.
    /// On an error `out` holds a partial document: the caller truncates it.
    fn query<F>(
        &mut self,
        run: &XmlRun<'_>,
        query: &XmlQuery<'_>,
        out: &mut Vec<u8>,
        hits: F,
    ) -> Result<()>
    where
        F: FnOnce(&mut HitSink<'_>) -> Result<()>;

    /// Appends the end of the document (the epilog), after the last query.
    fn finish(&mut self, out: &mut Vec<u8>) -> Result<()>;

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:757-763
    /// The message that NCBI posts (`ERR_POST` at the Error level, after the program's `Error: [name] `
    /// prefix) for a query whose report raised `exception` (`CaughtException`): the XML2 and JSON
    /// documents' `GetMsg()` and `what()`.
    fn caught_message(&self, exception: &crate::ncbi_exception::NcbiException) -> Vec<u8> {
        exception.what_after(exception.get_msg())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:361-369
    /// Appends what `PrintProlog` writes of the document, once: the start of the XML2 and JSON
    /// documents of outfmt 16 and 15 (nothing for BLAST XML, whose header comes with its first
    /// query, and for the side-file forms). NCBI writes it before it reads the first query
    /// batch, so a program's report calls this before the warnings of its first query, after
    /// its refusals of what NCBI could write (DECISIONS 20). `reported` is the run of the
    /// reports that the call writes after it (`None` for a call that writes none, a run without
    /// a query): its database counts that the reports would print where NCBI's sums overflow,
    /// which NCBI sums when it formats the first query's report (after its own errors of the
    /// search), are refused here, before the first byte (`common::database::check_sums_of`;
    /// BLASTX, which writes the start with the prolog, holds it, `blastx::blast_xml::HeldOutput`).
    fn start(&mut self, reported: Option<&XmlRun<'_>>, out: &mut Vec<u8>) -> Result<()>;

    /// Whether the document holds each query's report back until the next query's arrives, so
    /// that a run that stops (no epilog) never writes the report of the last query it formats:
    /// the XML2 and JSON documents of a `-subject` search (`Xml2Document`). BLAST XML writes each
    /// query's report at once.
    fn holds_back(&self) -> bool {
        false
    }

    /// Takes the place of the report of the last query of a run that stops, which the document
    /// holds back and never writes (`holds_back`): that report is not made, and `out` gets what
    /// the document writes when it arrives (the report held before it, unless the query's result
    /// joins it, `joins_held`: `XmlQuery::joins_held`).
    fn unwritten(&mut self, out: &mut Vec<u8>, joins_held: bool) -> Result<()> {
        let _ = (out, joins_held);
        anyhow::bail!("the document writes the report of every query")
    }

    /// What NCBI posted while the calls since the last `take_posts` wrote to their `out`
    /// (`DocumentPost`, at offsets of that buffer), in order, with the cut of the last query's
    /// loads and the exception that ends the run (`DocumentPosts`); the caller writes each post
    /// after the bytes before it (`report::query_warnings::write_hit_posts`), which hands the
    /// run's end back.
    fn take_posts(&mut self) -> DocumentPosts {
        DocumentPosts::default()
    }
}

/// How many of the `queries` of a report call have reports that the documents write, and so have
/// texts (titles, messages, the names of their hits' subjects) that NCBI's writer must be able
/// to write: all of them, except those of the last report when the run stops after this call
/// (`stops`: no epilog and no later report call follow it) and the documents hold that report back
/// (`holds_back`: `StructuredXml::holds_back`, `Xml2Run::holds_back`). That report starts at the
/// call's query `held_start` (the last query, or the first of the queries before it whose results
/// it joins, `XmlQuery::joins_held`; 0 where it starts in an earlier call).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2389
/// A `-subject` report is written when the next query's results arrive, or in the epilog
/// (`PrintEpilog`, blast_format.cpp:2206-2219), which a run that a later batch stops never
/// reaches: NCBI never makes the last formatted query's report, so it never fails on a title
/// there, and LOSAT refuses none of that report's texts.
pub fn written_queries(queries: usize, stops: bool, holds_back: bool, held_start: usize) -> usize {
    if stops && holds_back {
        held_start.min(queries.saturating_sub(1))
    } else {
        queries
    }
}

/// The hit-level calls of a structured writer behind `HitSink`.
pub trait HitWrite {
    fn hit(
        &mut self,
        hit: &XmlHit<'_>,
        hsps: &mut dyn FnMut(&mut HspSink<'_>) -> Result<()>,
    ) -> Result<()>;

    /// Where NCBI loads the subject of the next `Hit`, in the query's output
    /// (`HitSink::position`): where the `Hit` starts for a document that writes each `Hit` as it
    /// is made (outfmt 5), the start of the query's report for one that makes the whole report
    /// before it writes it (XML2 and JSON).
    fn position(&self) -> usize;
}

/// The HSP-level call of a structured writer behind `HspSink`.
pub trait HspWrite {
    fn hsp(&mut self, hsp: &XmlHsp<'_>) -> Result<()>;

    /// Whether the writer reads the coordinates, frames and rows of the HSPs of the `Hit` being
    /// written (`HspSink::reads_alignment`).
    fn reads_alignment(&self) -> bool {
        true
    }
}

/// The state of a document: the number of queries formatted so far (NCBI's
/// `SBlastXMLIncremental::m_IterationNum`), which numbers the iterations and tells whether the
/// header is still to be written.
pub struct BlastXmlWriter {
    iterations: u32,
    scratch: Scratch,
    posts: DocumentPosts,
}

impl BlastXmlWriter {
    /// A writer that has formatted `queries_before` queries (0 for a document's start; the
    /// count of the earlier calls when a program formats a document in several calls).
    pub fn new(queries_before: usize) -> Self {
        Self {
            iterations: queries_before as u32,
            scratch: Scratch::default(),
            posts: DocumentPosts::default(),
        }
    }
}

impl StructuredXml for BlastXmlWriter {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:745-760
    /// BLAST XML's message: NCBI formats one query per call (blast_format.cpp:708-740), so the
    /// label is always `Q(0/1` (without its `)`), then `what()`. The `Iteration` stays open after
    /// the `Hit`s written before the failed one (each written alone, `Xml5Hits`).
    fn caught_message(&self, exception: &crate::ncbi_exception::NcbiException) -> Vec<u8> {
        exception.what_after(b"Failed s_BlastXMLAddIteration Q(0/1")
    }

    /// Appends the `Iteration` of one query to `out` (before it the document's header, for the
    /// first query).
    fn query<F>(
        &mut self,
        run: &XmlRun<'_>,
        query: &XmlQuery<'_>,
        out: &mut Vec<u8>,
        hits: F,
    ) -> Result<()>
    where
        F: FnOnce(&mut HitSink<'_>) -> Result<()>,
    {
        use crate::report::ostream::{Cut, Site};
        if self.iterations == 0 {
            // The header is written outside the `try` of the queries (`Site::Report`).
            let header = out.len();
            if let Some(thrown) = write_header(out, run, query, &mut self.posts.posts)? {
                // No `Iteration` and no `Hit` of the query is made: the run ends after the posts.
                self.posts.ends_run = Some(Box::new(thrown.exception()));
                return Ok(());
            }
            // `serial_xml_start`, then the `<BlastOutput_iterations>` line: two insertions into
            // `out_stream`.
            self.posts.cuts.push(header, Cut::Site(Site::Report));
            self.posts.cuts.push(
                out.len() - b"\n<BlastOutput_iterations>".len(),
                Cut::Insertion,
            );
        }
        self.iterations += 1;
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:572-589
        // The whole `Iteration` without its hits is serialized first, cut after its length, and
        // the Hit objects are written one by one between the two parts. Its posts come before
        // it; an exception leaves the whole `Iteration` out (its number still counts) and no
        // hit is made.
        let start = out.len();
        let mut tail = std::mem::take(&mut self.scratch.tail);
        tail.clear();
        let mut enc = XmlEncoder::object(out);
        enc.begin_sequence("Iteration");
        int_member(&mut enc, "iter-num", 0, self.iterations as i32);
        string_member(&mut enc, "query-ID", 1, query.id);
        string_member(&mut enc, "query-def", 2, query_def(query.title));
        int_member(&mut enc, "query-len", 3, query.len);
        let head_end = enc.written();
        write_statistics(&mut enc, run, &query.stats);
        if let Some(message) = query.message {
            string_member(&mut enc, "message", 6, message);
        }
        enc.end_sequence();
        let iteration_posts = enc.take_posts();
        let thrown = enc.take_thrown();
        let error = enc.error();
        self.posts
            .posts
            .extend(replaced_posts(iteration_posts, run.app, |_| start));
        if let Some(error) = error {
            out.truncate(start);
            let Some(thrown) = thrown else {
                return Err(error.into());
            };
            self.posts
                .posts
                .push(iteration_post(start, run.app, &thrown));
            // No `Hit` of the query is made, so none of its subjects is loaded.
            self.posts.cut = true;
            self.scratch.tail = tail;
            return Ok(());
        }
        tail.extend_from_slice(&out[head_end..]);
        out.truncate(head_end);
        // The `Iteration`'s insertions, in the `try` (`Site::XmlIteration`): `serial_xml_start`,
        // `"\n<Iteration_hits>\n"`, each `Hit` (`Xml5Hits`), `"</Iteration_hits>"` and
        // `serial_xml_end`.
        // The `Iteration`'s regions, its parts and (most often) its `Hit`s, in one allocation.
        self.posts.cuts.reserve(32);
        self.posts.cuts.push(start, Cut::Site(Site::XmlIteration));
        self.posts.cuts.push(out.len(), Cut::Insertion);
        out.extend_from_slice(b"\n<Iteration_hits>\n");
        self.posts.cuts.push(out.len(), Cut::Insertion);
        let mut writer = Xml5Hits {
            run,
            masks: query.masks,
            out: &mut *out,
            scratch: &mut self.scratch,
            posts: &mut self.posts.posts,
            cuts: &mut self.posts.cuts,
            hits: 0,
        };
        let result = hits(&mut HitSink::new(&mut writer));
        match result {
            Ok(()) => {
                out.extend_from_slice(b"</Iteration_hits>");
                self.posts.cuts.push(out.len(), Cut::Insertion);
                out.extend_from_slice(&tail);
                self.posts.cuts.push(out.len(), Cut::Site(Site::Report));
            }
            // A `Hit` threw: the `Iteration` ends after the hits before it (no
            // `</Iteration_hits>`, no tail), no later `Hit` is made (nor its subject loaded), and
            // the run goes on.
            Err(error) if error.downcast_ref::<IterationThrew>().is_some() => {
                self.posts.cut = true;
                self.posts.cuts.push(out.len(), Cut::Site(Site::Report));
            }
            Err(error) => {
                // An exception that `s_BlastXMLAddIteration` catches (`CaughtException`) leaves
                // the `Iteration` open after the hits before it; no later `Hit` is made, nor its
                // subject loaded.
                if error.is::<super::CaughtException>() {
                    self.posts.cut = true;
                }
                self.scratch.tail = tail;
                return Err(error);
            }
        }
        self.scratch.tail = tail;
        Ok(())
    }

    fn take_posts(&mut self) -> DocumentPosts {
        std::mem::take(&mut self.posts)
    }

    /// The report-sum refusal of outfmt 5 (`XmlRun::check_sums`); BLAST XML writes nothing before
    /// its first query.
    fn start(&mut self, reported: Option<&XmlRun<'_>>, out: &mut Vec<u8>) -> Result<()> {
        let _ = out;
        reported.map_or(Ok(()), |run| run.check_sums(5))
    }

    /// Appends the end of the document if any query was formatted:
    /// `\n</BlastOutput_iterations>\n</BlastOutput>\n` and the epilog's end of line. With no
    /// query formatted the end of the document is empty and only the epilog's end of line
    /// remains.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2242-2243
    /// (`m_SerialXmlEnd` is set by the first formatted query, blastxml_format.cpp:755-763.)
    fn finish(&mut self, out: &mut Vec<u8>) -> Result<()> {
        use crate::report::ostream::{Cut, Site};
        // `m_SerialXmlEnd` and `endl` into `m_Outfile`: an insertion, the line feed and a flush, outside
        // the queries' `try`.
        self.posts.cuts.push(out.len(), Cut::Site(Site::Report));
        if self.iterations > 0 {
            out.extend_from_slice(b"\n</BlastOutput_iterations>\n</BlastOutput>\n");
            self.posts.cuts.push(out.len(), Cut::Insertion);
        }
        out.push(b'\n');
        self.posts.cuts.push(out.len(), Cut::Flush);
        Ok(())
    }
}

/// Where the hits of the query being written go: a structured writer's hit level
/// (`HitWrite`).
pub struct HitSink<'a> {
    inner: &'a mut dyn HitWrite,
}

impl<'a> HitSink<'a> {
    pub fn new(inner: &'a mut dyn HitWrite) -> Self {
        Self { inner }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:493-500
    /// Where NCBI loads the subject of the next `Hit`, in the query's output: what NCBI posts
    /// there (`report::query_warnings::HitPosts`) comes after the output before it. BLAST XML
    /// writes each `Hit` as NCBI makes it, and NCBI loads a `Hit`'s subject after the `Hit`s
    /// before it are written: the `Hit`'s start (the output's length). The XML2 and JSON
    /// documents make the whole report, loading every `Hit`'s subject, before they write any of
    /// it: the start of the report (`Xml2Document::report`).
    pub fn position(&self) -> usize {
        self.inner.position()
    }

    /// Writes one `Hit` (a subject with its HSPs, in NCBI's Seq-align order): `hsps` hands the
    /// HSPs one by one to the `HspSink`.
    pub fn hit<F>(&mut self, hit: &XmlHit<'_>, hsps: F) -> Result<()>
    where
        F: FnOnce(&mut HspSink<'_>) -> Result<()>,
    {
        let mut hsps = Some(hsps);
        self.inner.hit(hit, &mut |sink| {
            (hsps.take().expect("the writer asks for the HSPs once"))(sink)
        })
    }
}

/// Where the HSPs of the hit being written go: a structured writer's HSP level (`HspWrite`).
pub struct HspSink<'a> {
    inner: &'a mut dyn HspWrite,
}

impl<'a> HspSink<'a> {
    pub fn new(inner: &'a mut dyn HspWrite) -> Self {
        Self { inner }
    }

    /// Writes one `Hsp` (s_SeqAlignSetToXMLHsps, blastxml_format.cpp:107-339).
    pub fn hsp(&mut self, hsp: &XmlHsp<'_>) -> Result<()> {
        self.inner.hsp(hsp)
    }

    /// Whether the writer reads the alignment of each HSP of the `Hit` being written (its
    /// coordinates, frames and rows): not in outfmt 5's `Hit` of a subject whose lookup throws,
    /// nor in the `Hit` of a subject whose Bioseq NCBI's scope does not get (`XmlHit::unloaded`,
    /// outfmt 5, XML2 and JSON), whose HSPs NCBI writes with their Seq-align scores alone
    /// (`write_caught_hsp`). The program then
    /// hands each HSP as `XmlHsp::unaligned` (`HspSink::unaligned`), and reads nothing of the
    /// subject for it (its letters, its rows, the checks of its rows).
    pub fn reads_alignment(&self) -> bool {
        self.inner.reads_alignment()
    }

    /// Where the writer does not read the alignments of the `Hit` being written
    /// (`reads_alignment`): writes the HSPs `unaligned` (each `XmlHsp::unaligned`, made as they
    /// are written) and returns true, so that the program reads nothing of the subject for them;
    /// false, writing nothing, where the program writes each HSP with its alignment. The one
    /// answer of every program's outfmt 5 adapter (review 1 of line DBID pfoff-4, L3).
    pub fn unaligned<'h>(
        &mut self,
        unaligned: impl IntoIterator<Item = XmlHsp<'h>>,
    ) -> Result<bool> {
        if self.reads_alignment() {
            return Ok(false);
        }
        for hsp in unaligned {
            self.hsp(&hsp)?;
        }
        Ok(true)
    }
}

/// The hit level of outfmt 5: one `Hit` object per subject, written alone and joined into the
/// iteration.
struct Xml5Hits<'a> {
    run: &'a XmlRun<'a>,
    masks: &'a [QueryMask],
    out: &'a mut Vec<u8>,
    scratch: &'a mut Scratch,
    posts: &'a mut Vec<DocumentPost>,
    /// Each `Hit` is one insertion (`out_str` into `out_stream`).
    cuts: &'a mut crate::report::ostream::Cuts,
    hits: i32,
}

impl HitWrite for Xml5Hits<'_> {
    fn position(&self) -> usize {
        self.out.len()
    }

    fn hit(
        &mut self,
        hit: &XmlHit<'_>,
        hsps: &mut dyn FnMut(&mut HspSink<'_>) -> Result<()>,
    ) -> Result<()> {
        self.hits += 1;
        let scratch = &mut *self.scratch;
        let mut hit_buf = std::mem::take(&mut scratch.hit);
        hit_buf.clear();
        let result = (|| -> Result<()> {
            let mut enc = XmlEncoder::fragment(&mut hit_buf);
            enc.begin_document("Hit");
            enc.begin_sequence("Hit");
            int_member(&mut enc, "num", 0, self.hits);
            string_member(&mut enc, "id", 1, hit.id);
            string_member(&mut enc, "def", 2, query_def(hit.title));
            string_member(&mut enc, "accession", 3, hit.accession);
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:376-407
            // NCBI's catch writes 0: of a `-subject` record whose `GenerateDefline` threw
            // (`SubjectNames::of_unknown`), and of a database sequence whose lookup of its
            // Seq-align id throws (`SubjectNames::of_thrown`) or whose Bioseq the scope does not
            // get (`XmlHit::unloaded`), whose HSPs are the catch's too.
            let thrown = hit.lookup_throws.is_some() || hit.unloaded;
            let len = if hit.defline_thrown || thrown {
                0
            } else {
                hit.len
            };
            int_member(&mut enc, "len", 4, len);
            enc.member("hsps", 5);
            enc.begin_list(false);
            let mut writer = Xml5Hsps {
                run: self.run,
                masks: self.masks,
                hit: &mut enc,
                scratch,
                hsps: 0,
                thrown,
            };
            hsps(&mut HspSink::new(&mut writer))?;
            enc.end_list();
            enc.end_member();
            enc.end_sequence();
            enc.end_document();
            // NCBI writes each `Hit` into a string stream of its own (c++/src/algo/blast/format/blastxml_format.cpp:498-511)
            // and then to the output: its posts come before it, and an exception leaves it out
            // and ends the `Iteration` (`IterationThrew`).
            let at = self.out.len();
            if enc.has_posts() {
                self.posts
                    .extend(replaced_posts(enc.take_posts(), self.run.app, |_| at));
            }
            let Some(error) = enc.error() else {
                return Ok(());
            };
            match (error, enc.take_thrown()) {
                (_, Some(thrown)) => {
                    self.posts.push(iteration_post(at, self.run.app, &thrown));
                    Err(IterationThrew.into())
                }
                (error, None) => Err(error.into()),
            }
        })();
        if result.is_ok() {
            self.out.extend_from_slice(&hit_buf);
            self.cuts
                .push(self.out.len(), crate::report::ostream::Cut::Insertion);
        }
        scratch.hit = hit_buf;
        result
    }
}

/// The HSP level of outfmt 5. `thrown`: the `Hit`'s subject's lookup throws
/// (`XmlHit::lookup_throws`) or its Bioseq handle is empty (`XmlHit::unloaded`), so each HSP's
/// reads in NCBI's `try` throw too (`write_caught_hsp`).
struct Xml5Hsps<'a, 'o> {
    run: &'a XmlRun<'a>,
    masks: &'a [QueryMask],
    hit: &'a mut XmlEncoder<'o>,
    scratch: &'a mut Scratch,
    hsps: i32,
    thrown: bool,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:114-165
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:130-153
/// The `Hsp` number `number` that the `catch` of its subject's reads adds, in outfmt 5 (a subject
/// whose lookup throws, `XmlHit::lookup_throws`, or whose Bioseq NCBI's scope does not get,
/// `XmlHit::unloaded`) and, with `xml2`, in XML2 and JSON (the latter only, with their members'
/// tags and order): the number, the Seq-align's scores as `GetAlnScores` reads them (the bit score and
/// the e-value rounded to six digits as every HSP's), the Seq-align's `num_ident`
/// (`XmlHsp::alignment_identities`), the four coordinates 0 and empty rows; no frames, strands,
/// positives, gaps, length or midline.
pub(crate) fn write_caught_hsp<E: Encoder>(enc: &mut E, number: i32, hsp: &XmlHsp<'_>, xml2: bool) {
    let (identity, coordinates, rows) = if xml2 {
        (4, [9, 10, 13, 14], [19, 20])
    } else {
        (12, [4, 5, 6, 7], [17, 18])
    };
    enc.element();
    enc.begin_sequence("Hsp");
    int_member(enc, "num", 0, number);
    real_member(enc, "bit-score", 1, round6(hsp.bit_score));
    real_member(enc, "score", 2, f64::from(hsp.score));
    real_member(enc, "evalue", 3, round6(hsp.evalue));
    if xml2 {
        int_member(enc, "identity", identity, hsp.alignment_identities);
    }
    for (name, tag) in ["query-from", "query-to", "hit-from", "hit-to"]
        .into_iter()
        .zip(coordinates)
    {
        int_member(enc, name, tag, 0);
    }
    if !xml2 {
        int_member(enc, "identity", identity, hsp.alignment_identities);
    }
    string_member(enc, "qseq", rows[0], b"");
    string_member(enc, "hseq", rows[1], b"");
    enc.end_sequence();
}

/// The Dense-seg segments of an edit script as inclusive column ranges (`GetAlnStart` /
/// `GetAlnStop`); without a script, or with an empty one, one segment over the row (the whole
/// HSP as one substitution, as OF3's Seq-align writer does for an empty script,
/// `collect_seq_align_data`: `GapEditScript::Sub(q_end - q_offset)`).
fn fill_segments(segments: &mut Vec<(usize, usize)>, ops: Option<&[GapEditOp]>, width: usize) {
    segments.clear();
    match ops {
        Some(ops) if !ops.is_empty() => {
            let mut start = 0usize;
            for op in ops {
                let length = op.num() as usize;
                if length == 0 {
                    continue;
                }
                segments.push((start, start + length - 1));
                start += length;
            }
        }
        _ => {
            if width > 0 {
                segments.push((0, width - 1));
            }
        }
    }
}

/// What the XML and JSON writers count again from an HSP's rows (s_SeqAlignSetToXMLHsps,
/// blastxml_format.cpp:268-331, blastxml2_format.cpp:268-321): the masked query row
/// (`scratch.masked`: the regions of `masks` with the HSP's frame, in `mask_letter`), the
/// subject row as written (`HspCounts::hseq`), the midline (`scratch.midline`) and the three
/// counts. Both rows show their end gaps (`end_gaps_as_spaces`).
pub(crate) struct HspCounts {
    pub identity: i32,
    pub positive: i32,
    pub gaps: i32,
    /// Whether the subject row has end gaps, which `scratch.hseq` holds as spaces.
    subject_end_gaps: bool,
}

impl HspCounts {
    /// The subject row as written (`Hsp_hseq`, `"hseq"`): the HSP's row `srow` as it is, or,
    /// where it has end gaps (a Smith-Waterman alignment of `-use_sw_tback`), its copy with them
    /// as spaces (`scratch.hseq`, `count_hsp`).
    pub(crate) fn hseq<'r>(&self, srow: &'r [u8], scratch: &'r Scratch) -> &'r [u8] {
        if self.subject_end_gaps {
            &scratch.hseq
        } else {
            srow
        }
    }
}

pub(crate) fn count_hsp(
    hsp: &XmlHsp<'_>,
    positives: Option<&super::matrix::FormattingMatrix>,
    masks: &[QueryMask],
    mask_frame: i32,
    mask_letter: MaskLetter,
    scratch: &mut Scratch,
) -> HspCounts {
    fill_segments(&mut scratch.segments, hsp.ops, hsp.qrow.len());
    // The masked query row (only the regions of the HSP's own frame), then the recount over
    // the unmasked rows.
    mask_query_row(
        hsp.qrow,
        &scratch.segments,
        hsp.query_first,
        hsp.query_step,
        masks,
        mask_frame,
        mask_letter,
        &mut scratch.masked,
        &mut scratch.columns,
    );
    let (identity, positive) = recount(hsp.qrow, hsp.srow, positives, &mut scratch.midline);
    let ops = hsp.ops.unwrap_or_default();
    end_gaps_as_spaces(&mut scratch.masked, ops, AlignRow::Query);
    // The subject row is written as it is (`HspCounts::hseq`), copied only where it has end gaps.
    let subject_end_gaps =
        crate::report::serial::seqalign::row_end_gaps(ops, AlignRow::Subject).columns() > 0;
    if subject_end_gaps {
        scratch.hseq.clear();
        scratch.hseq.extend_from_slice(hsp.srow);
        end_gaps_as_spaces(&mut scratch.hseq, ops, AlignRow::Subject);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:238-241
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:225-230, 258-260
    // The gaps of the Seq-align's segments (also for a translated search, lines 259-260): the
    // columns of the edit script's gaps, not the `-` letters of the rows, which a stored gap
    // character of a protein subject also is. Without an edit script (an ungapped HSP: BLASTN's
    // `-ungapped` Dense-diag, which counts none, and TBLASTX's Std-seg, whose one Dense-seg segment
    // has no gap), none (remaining-gaps X1b-p).
    let gaps = hsp.ops.map_or(0, |ops| {
        ops.iter()
            .filter(|op| !op.is_sub())
            .map(|op| op.num() as usize)
            .sum::<usize>() as i32
    });
    HspCounts {
        identity,
        positive,
        gaps,
        subject_end_gaps,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:279-297
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:279-303
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:476-477
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:500-502
/// The row `row` (`AlignRow`) of an HSP whose edit script is `ops`, with its end gaps
/// (`serial::seqalign::row_end_gaps`: a Smith-Waterman alignment of `-use_sw_tback`) the end
/// character of `GetWholeAlnSeqString`, a space, as `Hsp_qseq` and `Hsp_hseq` show them. The
/// query row's masking and the midline read the row before (the masking counts the residues by
/// their `-`); the midline already has a space there. NCBI masks the row with its spaces and
/// skips only `-`, a mask that runs past the query's last letter reaching the row's last column;
/// the order is the same here because a Smith-Waterman alignment ends at its maximum, never with
/// a gap, and a mask starts at a letter, after the leading gaps.
fn end_gaps_as_spaces(written: &mut [u8], ops: &[GapEditOp], row: AlignRow) {
    let (lo, hi) = crate::report::serial::seqalign::row_end_gaps(ops, row).span(written.len());
    written[..lo].fill(b' ');
    written[hi..].fill(b' ');
}

impl HspWrite for Xml5Hsps<'_, '_> {
    fn reads_alignment(&self) -> bool {
        !self.thrown
    }

    fn hsp(&mut self, hsp: &XmlHsp<'_>) -> Result<()> {
        self.hsps += 1;
        if self.thrown {
            write_caught_hsp(&mut *self.hit, self.hsps, hsp, false);
            return Ok(());
        }
        let run = self.run;
        let scratch = &mut *self.scratch;
        let counts = count_hsp(
            hsp,
            run.positives,
            self.masks,
            hsp.query_frame,
            MaskLetter::Fixed(run.mask_char),
            scratch,
        );
        let enc = &mut *self.hit;
        enc.element();
        enc.begin_sequence("Hsp");
        int_member(enc, "num", 0, self.hsps);
        real_member(enc, "bit-score", 1, round6(hsp.bit_score));
        real_member(enc, "score", 2, f64::from(hsp.score));
        real_member(enc, "evalue", 3, round6(hsp.evalue));
        int_member(enc, "query-from", 4, hsp.query_from);
        int_member(enc, "query-to", 5, hsp.query_to);
        int_member(enc, "hit-from", 6, hsp.hit_from);
        int_member(enc, "hit-to", 7, hsp.hit_to);
        int_member(enc, "query-frame", 10, hsp.query_frame);
        int_member(enc, "hit-frame", 11, hsp.hit_frame);
        int_member(enc, "identity", 12, counts.identity);
        int_member(enc, "positive", 13, counts.positive);
        int_member(enc, "gaps", 14, counts.gaps);
        int_member(enc, "align-len", 15, hsp.qrow.len() as i32);
        string_member(enc, "qseq", 17, &scratch.masked);
        string_member(enc, "hseq", 18, counts.hseq(hsp.srow, scratch));
        string_member(enc, "midline", 19, &scratch.midline);
        enc.end_sequence();
        Ok(())
    }
}

#[cfg(test)]
#[path = "writer_segment_tests.rs"]
mod segment_tests;
