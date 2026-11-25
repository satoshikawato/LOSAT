//! Warnings that NCBI writes to standard error for a search, whatever the output format.

use crate::blastinput::fasta_reader::InputRecord;
use crate::report::archive::ArchiveHandle;
use crate::report::diag_log::Posts;
use anyhow::Result;
use std::io::{self, Write};

/// The warnings that NCBI posts between the reports of a search: the title warnings of a
/// query batch when it is read (before the report of its first query) and the warning of
/// an invalid query when its report starts (`PrintOneResultSet`, before the query's
/// preamble). NCBI posts them on `cerr`, which the C++ library ties to `cout`, so the
/// report written so far reaches standard output before each of them; a failed write of
/// the report stops NCBI there (its stream throws, `blast_format.cpp:118-119`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:4086
/// NCBI reference (598d8ae6): c++/include/corelib/ncbistre.hpp:543-544
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
pub struct QueryWarnings<'a> {
    /// The posts to write before the report of each query, by query index (`of_posts`).
    before: &'a [Posts],
    /// The BLAST archive's handle, which saves each whole post as one message (`Posts::write`).
    archive: Option<&'a ArchiveHandle>,
    /// Standard error (the CLI) or the caller's diagnostics.
    pub sink: &'a mut dyn Write,
    /// The texts of the subjects that the reports read one by one (`with_subject_reads`), and
    /// how many of them are written.
    reads: Option<ReadCursor<'a>>,
    /// The posts of each query whose later lines (`Posts::later`) this format's reports write,
    /// and the format's receiver of them (`with_later`).
    later: Option<(&'a [Posts], &'a mut Vec<u8>)>,
    /// Where a `-html` report's own posts go apart from `sink` (`with_html_posts`), `None`: to
    /// `sink`.
    html_posts: Option<&'a mut dyn Write>,
    /// The report's place among its search's formats (`for_report`), which keeps the first call
    /// of each report's id links (`formatter_registry::InitConfig::call_posts`).
    report: usize,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1587-1589
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1947-1955
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:374-375
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1465-1466
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

/// The reads of one query's report (`SubjectReads`), in the order of the reads.
pub type QueryReads = Vec<SubjectRead>;

/// One read of a query's report: the subject as the report's writer names it (its index in the
/// hit lists), the text that its load writes in the program's model, and, in a run whose lookups
/// can write a text, what the load looks up (`ReadLoad`; then every read has an entry, with or
/// without a text of its own).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubjectRead {
    pub subject: u32,
    pub text: Vec<u8>,
    pub load: Option<ReadLoad>,
}

impl SubjectRead {
    /// A read of `subject` whose text is `text`, without the facts of its load (the tests of the
    /// writers).
    #[cfg(test)]
    pub fn of_text(subject: u32, text: &[u8]) -> Self {
        Self {
            subject,
            text: text.to_vec(),
            load: None,
        }
    }

    /// The text of this read where it happens: the texts of the loads left out before that this
    /// load makes ([`ReadLoad::placed`]) around the model's own (also a load of a batch's prefetch,
    /// `diag_log::Posts::loads`).
    pub fn placed(&self) -> Vec<u8> {
        self.placed_in(None)
    }

    /// [`SubjectRead::placed`] with the texts left out in `ledger` where it is given instead of
    /// the run's ([`LoadTexts::left_before`], the model of the writers that walks the writers'
    /// own steps in a ledger of its own).
    pub fn placed_in(&self, ledger: Option<&std::sync::Arc<LoadTexts>>) -> Vec<u8> {
        match (&self.load, ledger) {
            (Some(load), Some(ledger)) => load.in_ledger(ledger).placed(&self.text),
            (Some(load), None) => load.placed(&self.text),
            (None, _) => self.text.clone(),
        }
    }

    /// [`ReadLoad::place`] of this read (a load of a batch's prefetch that stops,
    /// `common::database::prefetch_stop`) with the texts left out `left` instead of the run's
    /// ledger: those that the reports and prefetches before it leave where they are written
    /// ([`LoadTexts::left_before`]), for a program that models the batch before it writes them.
    pub fn place_left(&mut self, left: &std::sync::Arc<LoadTexts>) {
        if let Some(load) = &mut self.load {
            let mut moved = load.in_ledger(left);
            self.text = moved.place(std::mem::take(&mut self.text));
            (load.lead, load.opening) = (moved.lead, moved.opening);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1440-1466 (cited at
// `TaxdbRun::load`): a load asks the taxonomy database for a taxid that its volume's cache does
// not hold, and the run's first ask opens it (the static `CTaxDBFileInfo` of `GetTaxNames`).
/// What the load of a read's subject looks up, as far as it decides the texts: NCBI's lookups write
/// the texts of opening the taxonomy database at the run's first ask, and the offset error at a
/// volume's first lookup of the taxid whose record begins past the end of `taxdb.btd`
/// (`TaxDb::offset_error_taxid`; on a database with a volume whose cache can ask again,
/// `TaxdbRun::cache_collision`, at each lookup that misses the slot of that taxid), nothing else.
/// The model of the reads gives each text to the load that makes the ask, assuming every load it
/// models happens; a load that NCBI does not make (a BLAST XML `Iteration` or `Hit` that the stream
/// leaves out, `ReadCursor::cut`) leaves its texts to NCBI's next load that makes the same ask
/// (`LoadTexts`), and is LOSAT's refusal where it takes that slot on such a database
/// (`cut_refusal`: the model followed the slot through it).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadLoad {
    /// The run's texts of the loads left out (shared by every read of the run).
    pub texts: std::sync::Arc<LoadTexts>,
    /// The volume of the loaded blob.
    pub volume: u32,
    /// Whether the load looks a positive taxid up (so it makes the run's first ask, where no load
    /// before it asked).
    pub looks_up: bool,
    /// Whether the load looks up the taxid of the offset error.
    pub offset: bool,
    /// The bytes of the model's text that are NCBI's warning of a UTF8String read for a
    /// VisibleString (its start): the lookup of the load's Seq-align id reads such a header first
    /// (`ParsedSubject::scope_utf8`, before the blob's taxonomy), and the run posts it here.
    pub lead: usize,
    /// The bytes of the model's text after `lead` that are the texts of opening.
    pub opening: usize,
    /// Whether the lookup of the load's Seq-align id reads such a header: the load takes the
    /// warning that a load NCBI does not make left (`left_out`).
    pub scope_utf8: bool,
    /// LOSAT's refusal where NCBI does not make this load (`left_out`), whose lookups take the
    /// slot of the offset error's taxid on a database whose volume holds that taxid and another of
    /// its slot (`common::database::TaxdbRun::cache_collision`): the model of the later loads
    /// followed the slot through this one. `None` for any other load.
    pub cut_refusal: Option<std::sync::Arc<str>>,
}

/// LOSAT's refusal of a load that NCBI's BLAST XML stream leaves out where LOSAT's model of the
/// taxonomy cache's slot followed the slot through it ([`ReadLoad::cut_refusal`]): no error of
/// NCBI's output stream. It travels as the `io::Error` of the report's writes
/// ([`CutRefusal::error`]), and a program that turns its report's write errors into NCBI's takes
/// it out as LOSAT's refusal ([`CutRefusal::of`]), which drops what the run holds.
#[derive(Debug)]
pub(crate) struct CutRefusal(std::sync::Arc<str>);

impl CutRefusal {
    /// The refusal `text` as the error of the report's writes.
    fn error(text: &std::sync::Arc<str>) -> io::Error {
        io::Error::new(
            io::ErrorKind::Unsupported,
            CutRefusal(std::sync::Arc::clone(text)),
        )
    }

    /// The refusal's text, where `error` is a [`CutRefusal`] of the report's writes.
    pub(crate) fn of(error: &io::Error) -> Option<&str> {
        error
            .get_ref()
            .and_then(|inner| inner.downcast_ref::<CutRefusal>())
            .map(|refusal| &*refusal.0)
    }
}

impl std::fmt::Display for CutRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for CutRefusal {}

impl ReadLoad {
    /// The load of a subject whose Seq-align id's lookup reads a header with a UTF8String where a
    /// VisibleString is expected and whose blob's taxonomy writes no text (no `taxdb.bti`, or no
    /// blob), with the run's ledger `texts`.
    pub fn of_scope_utf8(texts: std::sync::Arc<LoadTexts>) -> Self {
        ReadLoad {
            texts,
            volume: u32::MAX,
            looks_up: false,
            offset: false,
            lead: 0,
            opening: 0,
            scope_utf8: true,
            cut_refusal: None,
        }
    }

    /// The text that this load writes where it happens, `own` being the model's: after its
    /// warning of a UTF8String header (`lead`), or where it has none and its lookup reads such a
    /// header, the one that a load left out (`left_out`), the texts of opening that a load left
    /// out placed, where this load looks a taxid up, and after it the offset error that a load of
    /// the same volume left out, where this load looks that taxid up (the texts of a load are its
    /// lookup's warning, then those of opening, then those of its lookups). They are NCBI's texts
    /// of this load: the asks that the left-out load made in the model are this load's first.
    pub fn placed(&self, own: &[u8]) -> Vec<u8> {
        self.with_left_out(own.to_vec()).0
    }

    /// `placed` for the model's read of this load (`common::database::subject_read_texts`, a
    /// batch modelled after a report that left loads out): its text, whose warning and texts of
    /// opening are then those placed before its own.
    pub fn place(&mut self, own: Vec<u8>) -> Vec<u8> {
        let (text, lead, opening) = self.with_left_out(own);
        (self.lead, self.opening) = (lead, opening);
        text
    }

    /// `place` of the model's text `own` after the warning `lead` that this load's lookup posts
    /// (empty where it posts none).
    pub fn place_after(&mut self, lead: Vec<u8>, own: Vec<u8>) -> Vec<u8> {
        self.lead = lead.len();
        let mut text = lead;
        text.extend(own);
        self.place(text)
    }

    /// The text of `placed`, the length of its warning and that of its texts of opening.
    fn with_left_out(&self, own: Vec<u8>) -> (Vec<u8>, usize, usize) {
        let mut pending = self.texts.pending.lock().expect("the texts of the loads");
        if pending.is_empty() {
            return (own, self.lead, self.opening);
        }
        let (lead, own) = own.split_at(self.lead.min(own.len()));
        let mut text = lead.to_vec();
        if lead.is_empty() && self.scope_utf8 {
            text.append(&mut pending.utf8);
        }
        let lead = text.len();
        if self.looks_up {
            text.append(&mut pending.opening);
        }
        let opening = text.len() - lead + self.opening;
        text.extend_from_slice(own);
        if self.offset {
            if let Some(at) = pending
                .offsets
                .iter()
                .position(|(volume, _)| *volume == self.volume)
            {
                text.extend(pending.offsets.remove(at).1);
            }
        }
        (text, lead, opening)
    }

    /// This load with the texts left out in `texts` instead of the run's ledger
    /// ([`LoadTexts::left_before`]).
    fn in_ledger(&self, texts: &std::sync::Arc<LoadTexts>) -> ReadLoad {
        ReadLoad {
            texts: std::sync::Arc::clone(texts),
            ..self.clone()
        }
    }

    /// NCBI does not make this load, whose text in the model is `text`: its warning of a
    /// UTF8String header goes to the next load whose lookup reads such a header, its texts of
    /// opening and its offset error to the next load that makes those asks (`placed`). `Err`:
    /// LOSAT's refusal where the model followed the taxonomy cache's slot through this load
    /// (`cut_refusal`).
    pub fn left_out(&self, text: &[u8]) -> io::Result<()> {
        if let Some(refusal) = &self.cut_refusal {
            return Err(CutRefusal::error(refusal));
        }
        let (lead, text) = text.split_at(self.lead.min(text.len()));
        let (opening, offset) = text.split_at(self.opening.min(text.len()));
        let mut pending = self.texts.pending.lock().expect("the texts of the loads");
        pending.utf8.extend_from_slice(lead);
        if !opening.is_empty() {
            pending.opening.extend_from_slice(opening);
        }
        if !offset.is_empty() {
            match pending
                .offsets
                .iter_mut()
                .find(|(volume, _)| *volume == self.volume)
            {
                Some((_, texts)) => texts.extend_from_slice(offset),
                None => pending.offsets.push((self.volume, offset.to_vec())),
            }
        }
        Ok(())
    }
}

/// The texts of the loads of a run that NCBI does not make (`ReadLoad::left_out`), until the load
/// that makes their asks (`ReadLoad::placed`): a load of a later read of the report, of a later
/// batch's reads, or a later batch's prefetch (`common::database::subject_read_texts`). A text
/// that no later load asks for is never written. One per run (`TaxdbRun`).
#[derive(Debug, Default)]
pub struct LoadTexts {
    pending: std::sync::Mutex<PendingTexts>,
}

/// The run's one ledger: two are the same where they are one (the loads of one run share it).
impl PartialEq for LoadTexts {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self, other)
    }
}

impl Eq for LoadTexts {}

impl LoadTexts {
    /// No text left out (`TaxdbRun::restart`: the batches are searched again).
    pub fn clear(&self) {
        *self.pending.lock().expect("the texts of the loads") = PendingTexts::default();
    }

    /// Whether no text is left out.
    pub fn is_empty(&self) -> bool {
        self.pending
            .lock()
            .expect("the texts of the loads")
            .is_empty()
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288 (cited at
    // `query_posts`: a batch's prefetch, then the result sets of its queries)
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:372-375 (cited at
    // `ReadCursor::cut`: a report's `Hit` loads its subject)
    /// The texts left out, in a ledger of their own, where the writers of a run that left nothing
    /// before have written the reports of the queries below `query` and the posts before each:
    /// the writers' own steps, run here before any of it is written, in that ledger
    /// (`ReadCursor::in_ledger`). By query, the step before its report (`before_report`, the
    /// writers' `QueryWarnings::before_query`: the reads of the reports before that they did not
    /// place, then the posts written before its report, `prefetched`, a batch's prefetch at its
    /// first query, `Posts::loads`), then its report: one whose outfmt 5 `Iteration` NCBI's XML
    /// stream cut before the first `Hit` (`cut`) makes none of its reads (`ReadCursor::cut`, as
    /// the writers' `write_posts_to` at that cut); another makes them all, at its `Hit`s or after
    /// it (`ReadCursor::subject_read`, `ReadCursor::before_query`; the ledger sees them in the same
    /// order either way). Last, the reads of the report before `query` that it did not place. For
    /// a program that models a later batch's prefetch before it writes these reports
    /// (`common::database::XmlReportsBefore`): what the run's ledger holds at that prefetch where
    /// they are written (review 6 of line DBID loadord, L2: one walk). A cut at a later `Hit` is
    /// the writer's alone (the model writes no report); a database title cannot make it (review 6).
    pub fn left_before<'p>(
        reads: &[QueryReads],
        cut: &dyn Fn(usize) -> bool,
        prefetched: &dyn Fn(usize) -> Option<&'p Posts>,
        query: usize,
    ) -> io::Result<std::sync::Arc<LoadTexts>> {
        let left = std::sync::Arc::new(LoadTexts::default());
        let mut cursor = ReadCursor::in_ledger(reads, &left);
        // Nothing is written: the model keeps the ledger only.
        let (mut report, mut sink) = (io::sink(), io::sink());
        for before in 0..query {
            // The posts' texts are made in the ledger; the model writes none of them.
            let _placed = before_report(
                Some(&mut cursor),
                prefetched(before),
                before,
                &mut report,
                &mut sink,
            )
            .expect("a model writes to no stream");
            if cut(before) {
                // The only error of a model: LOSAT's refusal of a load left out (`left_out`).
                cursor.cut(before, &mut report, &mut sink)?;
            }
        }
        cursor
            .before_query(query, &mut report, &mut sink)
            .expect("a model writes to no stream");
        Ok(left)
    }
}

#[derive(Debug, Default)]
struct PendingTexts {
    /// NCBI's warning of a UTF8String header that a left-out load's lookup read (`ReadLoad::lead`).
    utf8: Vec<u8>,
    /// The texts of opening the taxonomy database (empty: none left out).
    opening: Vec<u8>,
    /// By volume, its offset error.
    offsets: Vec<(u32, Vec<u8>)>,
}

impl PendingTexts {
    fn is_empty(&self) -> bool {
        self.utf8.is_empty() && self.opening.is_empty() && self.offsets.is_empty()
    }
}

impl SubjectReads {
    /// Sets the reads of query `query`.
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
                .filter_map(|read| {
                    let rank = ordinals.binary_search(&read.subject).ok()?;
                    Some(SubjectRead {
                        subject: rank as u32,
                        ..read.clone()
                    })
                })
                .collect()
        })
        .collect()
}

/// The reads of query `query` in `reads`, in order.
fn reads_of(reads: &[QueryReads], query: usize) -> &[SubjectRead] {
    reads.get(query).map_or(&[], Vec::as_slice)
}

/// A report writer's call where NCBI reads a subject (`QueryWarnings::subject_read`): the
/// subject's index and the writer, flushed before a text is written or, for a writer that
/// formats into a buffer, marked at the buffer's length.
pub type SubjectReadHook<'r, W> = &'r mut dyn FnMut(u32, &mut W) -> io::Result<()>;

/// A place of a pairwise report (outfmt 0 to 4) where NCBI writes to the run's diagnostics between
/// the report's bytes, after the report written before it (`cerr`, tied to `cout`): the load of a
/// subject (`Load`: its read's text, `QueryWarnings::subject_read`) or, with `-html`, a read of
/// NCBI's formatter registries (`Html`: `report::html_ids::HtmlPost`, whose text the run's
/// `-html` links give, `HtmlLinks::posts_of`). The writers call their `ReportHook` at each, in the
/// report's order; one that formats into a buffer keeps them with their offsets
/// (`x_pairwise_par::format_piece`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReportPlace {
    Load(u32),
    Html(crate::report::html_ids::HtmlPost),
}

/// The pairwise writers' call at each `ReportPlace` of a query's report, with the report written so
/// far (`x_pairwise_par::at_place`, or a writer's own marks).
pub type ReportHook<'r, W> = &'r mut dyn FnMut(ReportPlace, &mut W) -> io::Result<()>;

/// A [`ReportHook`] that writes nothing: of a writer without `-html` that reads no subject
/// through it (psiblast's headings), and of the tests of the writers.
pub fn no_places<W>(_place: ReportPlace, _writer: &mut W) -> io::Result<()> {
    Ok(())
}

/// The writing of the texts of `SubjectReads` where the reports read the subjects: every text
/// of the queries below `query` is written, and the first `next` of query `query`'s. A read
/// writes its text where it happens (`SubjectRead::placed`), with the texts that the loads NCBI
/// does not make left to it (`cut`).
pub struct ReadCursor<'a> {
    reads: &'a [QueryReads],
    query: usize,
    next: usize,
    /// The run ended at an exception of NCBI's XML stream (`end_run`): no text is written after.
    ended: bool,
    /// The ledger of the texts left out where it is not the run's (`in_ledger`: the model of the
    /// writers, `LoadTexts::left_before`); `None` for the writers.
    ledger: Option<std::sync::Arc<LoadTexts>>,
}

impl<'a> ReadCursor<'a> {
    /// The cursor of `reads` (by query index), or `None` where no query has a read.
    pub fn new(reads: &'a [QueryReads]) -> Option<Self> {
        reads
            .iter()
            .any(|query| !query.is_empty())
            .then_some(ReadCursor {
                reads,
                query: 0,
                next: 0,
                ended: false,
                ledger: None,
            })
    }

    /// The cursor of `reads` whose loads take and leave the texts left out in `ledger` instead of
    /// the run's: the writers' steps run by the model of the writers before they write
    /// (`LoadTexts::left_before`).
    fn in_ledger(reads: &'a [QueryReads], ledger: &std::sync::Arc<LoadTexts>) -> Self {
        ReadCursor {
            reads,
            query: 0,
            next: 0,
            ended: false,
            ledger: Some(std::sync::Arc::clone(ledger)),
        }
    }

    /// The text of `read` where it happens (`SubjectRead::placed`, in the cursor's ledger).
    fn placed(&self, read: &SubjectRead) -> Vec<u8> {
        read.placed_in(self.ledger.as_ref())
    }

    /// Whether the report of query `query` reads a subject that may write a text.
    pub fn has_texts(&self, query: usize) -> bool {
        !self.ended && !reads_of(self.reads, query).is_empty()
    }

    /// The report of query `query` reads subject `subject` here (a tabular row that reads it,
    /// its description line, its alignment's heading, its BLAST XML `Hit`): the texts of the
    /// query's reads up to its first are written to `sink`, after the report written so far
    /// (`report`) is flushed, unless they are written already. A subject that the query's report
    /// read before, or that is not among its reads, writes nothing.
    pub fn subject_read(
        &mut self,
        query: usize,
        subject: u32,
        report: &mut dyn Write,
        sink: &mut dyn Write,
    ) -> io::Result<()> {
        if query < self.query || self.ended {
            return Ok(());
        }
        let from = if query == self.query { self.next } else { 0 };
        let reads = reads_of(self.reads, query);
        let Some(found) = reads
            .get(from..)
            .and_then(|rest| rest.iter().position(|read| read.subject == subject))
        else {
            return Ok(());
        };
        let mut text = self.unwritten_before(query);
        for read in &reads[from..=from + found] {
            text.extend(self.placed(read));
        }
        self.next = from + found + 1;
        write_after(report, sink, &text)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:372-375
    /// The report of query `query` made no read after the ones written (BLAST XML: its `Iteration`
    /// threw, or one of its `Hit`s did, `blast_xml::DocumentPosts::cut`): NCBI loads a subject
    /// when it makes its `Hit` (`s_SeqAlignToXMLHit`), so the query's later reads do not happen.
    /// Their texts go to NCBI's next loads that make the same asks (`ReadLoad::left_out`); one
    /// that no later load asks for is never written. The texts of the queries before it that
    /// their reports did not place are written first, to `sink` after `report` is flushed
    /// (`before_query`).
    pub fn cut(
        &mut self,
        query: usize,
        report: &mut dyn Write,
        sink: &mut dyn Write,
    ) -> io::Result<()> {
        if query < self.query || self.ended {
            return Ok(());
        }
        self.before_query(query, report, sink)?;
        let reads = reads_of(self.reads, query);
        for read in &reads[self.next.min(reads.len())..] {
            match (&read.load, &self.ledger) {
                (Some(load), Some(ledger)) => load.in_ledger(ledger).left_out(&read.text)?,
                (Some(load), None) => load.left_out(&read.text)?,
                (None, _) => {}
            }
        }
        self.next = reads.len();
        Ok(())
    }

    /// The run ends at an exception of NCBI's XML stream (`blast_xml::DocumentPosts::ends_run`):
    /// NCBI loads nothing after it, so no text not written yet is written.
    pub fn end_run(&mut self) {
        self.ended = true;
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
        if self.ended {
            return Ok(());
        }
        let text = self.unwritten_before(query);
        write_after(report, sink, &text)
    }

    /// The texts not written yet of every query (a report that ends without placing them).
    fn unwritten(&mut self) -> Vec<u8> {
        if self.ended {
            return Vec::new();
        }
        let queries = self.reads.len().max(self.query + 1);
        self.unwritten_before(queries)
    }

    /// The texts not written yet of the queries below `query`, which count as written (their
    /// loads as made) after this; the cursor moves to the start of `query`'s.
    fn unwritten_before(&mut self, query: usize) -> Vec<u8> {
        if query <= self.query {
            return Vec::new();
        }
        let mut text = Vec::new();
        for read in reads_of(self.reads, self.query).iter().skip(self.next) {
            text.extend(self.placed(read));
        }
        for later in self.query + 1..query {
            for read in reads_of(self.reads, later) {
                text.extend(self.placed(read));
            }
        }
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
        ..HitPosts::default()
    };
    match warnings {
        Some(warnings) => {
            let QueryWarnings { reads, sink, .. } = warnings;
            write_reads(report, &posts, query, writer, &mut **sink, reads.as_mut())
        }
        // Marks alone post nothing of their own.
        None => writer.write_all(report),
    }
}

/// The places in a query's report, formatted into a buffer, where NCBI loads a subject and posts
/// texts: a structured report's `Hit` (BLAST XML, XML2, JSON; `HitSink::position`) or a tabular
/// row's or pairwise line's read (`write_marked`). At each, the text of the subject's read in the
/// program's model (`ReadCursor::subject_read`) and the texts that the report's writer makes there
/// (the taxonomy database's messages of the XML2 and JSON documents' `Hit`s, `TaxdbPrinter`).
/// NCBI posts them on `cerr`, which flushes the report written before them (`QueryWarnings`);
/// `write` writes them so. Dropping it loses them and the exception that ends the run, so it is
/// `must_use`.
#[derive(Debug, Default)]
#[must_use = "a report's posts and the exception that ends the run go through write_hit_posts"]
pub struct HitPosts {
    /// Whether every `Hit` is marked (the model's reads of the query write texts); else only one
    /// whose writer makes a text.
    every: bool,
    /// By marked `Hit`: its place in the query's output and its subject.
    marks: Vec<(usize, u32)>,
    /// By marked `Hit`: the end of its texts in `texts`.
    ends: Vec<usize>,
    texts: Vec<u8>,
    /// What the document's writer posted while it wrote the report
    /// (`StructuredXml::take_posts`), at offsets of the report, in order.
    document: Vec<crate::report::blast_xml::writer::DocumentPost>,
    /// The report made no read after its last marked `Hit` (`DocumentPosts::cut`).
    cut: bool,
    /// The program's post of an exception that the report raised and NCBI's report function
    /// caught (`blast_xml::writer::caught_post`), after every other post of the report.
    caught: Vec<u8>,
    /// The exception that ends the run after these posts (`DocumentPosts::ends_run`).
    ends_run: Option<Box<crate::ncbi_exception::NcbiException>>,
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

    /// The `Hit` `hit` of subject `subject` loads it at `at` (`HitSink::position`); `post` writes
    /// the texts that the writer makes there. A `Hit` whose subject's lookup throws
    /// (`XmlHit::lookup_throws`) loads nothing: no read and no text of it (`post` is not called).
    pub fn hit(
        &mut self,
        at: usize,
        hit: &crate::report::blast_xml::XmlHit<'_>,
        subject: u32,
        post: impl FnOnce(&mut Vec<u8>) -> io::Result<()>,
    ) -> io::Result<()> {
        if hit.lookup_throws.is_some() {
            return Ok(());
        }
        self.mark(at, subject, post)
    }

    /// `hit` of a `Hit` that loads its subject.
    fn mark(
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

    /// Adds what the document's writer posted while it wrote the report (`DocumentPost`, at
    /// offsets of the report): each comes after the report's bytes before it, and after a `Hit`'s
    /// texts at the same place (NCBI loads a subject before it writes the `Hit`); with the cut of
    /// the report's loads and the exception that ends the run after them (`DocumentPosts`).
    pub fn document(&mut self, posts: crate::report::blast_xml::DocumentPosts) {
        let crate::report::blast_xml::DocumentPosts {
            posts,
            cut,
            ends_run,
        } = posts;
        self.document.extend(posts);
        self.cut |= cut;
        if ends_run.is_some() {
            self.ends_run = ends_run;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:757-760
    /// The report raised an exception that NCBI's report function catches: its post `post`
    /// (`blast_xml::writer::caught_post`) comes after the texts of the loads that the report made
    /// before it and the document's posts (`write_posts_to`).
    pub fn caught(&mut self, post: Vec<u8>) {
        self.caught = post;
    }

    /// Whether nothing of these posts goes to the diagnostics (no text, no document post).
    fn posts_nothing(&self) -> bool {
        self.texts.is_empty() && self.document.is_empty() && self.caught.is_empty()
    }

    /// Writes `report`, the query's output, to `writer`, and at each marked `Hit`, after the
    /// bytes before it: `read` of its subject (the model's text, to `sink`; `writer` is the
    /// report) and the texts the writer made there, to `sink` after `writer` is flushed; the
    /// document's posts (`document`) likewise at their places.
    pub fn write(
        &self,
        report: &[u8],
        writer: &mut dyn Write,
        sink: &mut dyn Write,
        mut read: impl FnMut(u32, &mut dyn Write, &mut dyn Write) -> io::Result<()>,
    ) -> io::Result<()> {
        let (mut index, mut from) = (0, 0);
        let mut document = self.document.iter().peekable();
        let mut written = 0;
        for &(at, subject) in &self.marks {
            while let Some(post) = document.next_if(|post| post.at < at) {
                writer.write_all(&report[written..post.at])?;
                written = post.at;
                write_after(writer, &mut *sink, &post.text)?;
            }
            writer.write_all(&report[written..at])?;
            written = at;
            read(subject, writer, &mut *sink)?;
            let end = self.ends[index];
            index += 1;
            let texts = &self.texts[from..end];
            from = end;
            write_after(writer, &mut *sink, texts)?;
        }
        for post in document {
            writer.write_all(&report[written..post.at])?;
            written = post.at;
            write_after(writer, &mut *sink, &post.text)?;
        }
        writer.write_all(&report[written..])
    }
}

/// Writes `report`, a part of a structured document that holds no `Hit` of its own (the end of the
/// document, a held report written there), with what the document's writer posted while it wrote
/// it (`StructuredXml::take_posts`), as `write_hit_posts`.
pub fn write_document_posts(
    report: &[u8],
    posts: crate::report::blast_xml::DocumentPosts,
    writer: &mut dyn Write,
    warnings: Option<&mut QueryWarnings<'_>>,
) -> Result<()> {
    let mut hits = HitPosts::default();
    hits.document(posts);
    write_hit_posts(report, &hits, 0, writer, warnings)
}

/// Writes query `query`'s structured report `report` with the posts of its `Hit`s
/// (`HitPosts::write`): the model's reads through the warnings' cursor
/// (`QueryWarnings::subject_read`) and the writer's texts, both to the warnings' sink (the run's
/// diagnostics; `write_posts_to`). Without warnings, the report alone: a post there would be lost,
/// which is LOSAT's internal error.
pub fn write_hit_posts(
    report: &[u8],
    hits: &HitPosts,
    query: usize,
    writer: &mut dyn Write,
    warnings: Option<&mut QueryWarnings<'_>>,
) -> Result<()> {
    match warnings {
        Some(warnings) => {
            let QueryWarnings { reads, sink, .. } = warnings;
            write_posts_to(report, hits, query, writer, &mut **sink, reads.as_mut())
        }
        None => {
            anyhow::ensure!(
                hits.posts_nothing(),
                "LOSAT internal error: a structured report posts on the diagnostics, but its format writes no warnings"
            );
            write_posts_to(report, hits, query, writer, &mut std::io::sink(), None)
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:775-777
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
/// The one writer of a structured report's posts (`write_hit_posts`; BLASTX, whose reports keep
/// their own cursor, calls it with that): `report` to `writer` with the posts of its `Hit`s and of
/// its document (`HitPosts::write`) to `sink`, the reads' texts through `reads`. A report that made
/// no read after its last `Hit` (`DocumentPosts::cut`) moves the texts of the later reads
/// (`ReadCursor::cut`); the post of an exception that the report function caught
/// (`HitPosts::caught`) comes after them; an exception of NCBI's XML stream that ends the run is
/// returned after the posts (`blast_xml::writer::run_ending_exception`), and no text of a read is
/// written after it.
pub fn write_posts_to(
    report: &[u8],
    hits: &HitPosts,
    query: usize,
    writer: &mut dyn Write,
    sink: &mut dyn Write,
    mut reads: Option<&mut ReadCursor<'_>>,
) -> Result<()> {
    write_reads(
        report,
        hits,
        query,
        writer,
        &mut *sink,
        reads.as_deref_mut(),
    )?;
    if let Some(cursor) = reads {
        if hits.ends_run.is_some() {
            cursor.end_run();
        } else if hits.cut {
            cursor.cut(query, writer, &mut *sink)?;
        }
    }
    write_after(writer, sink, &hits.caught)?;
    match &hits.ends_run {
        Some(thrown) => Err(crate::report::blast_xml::writer::run_ending_exception(
            thrown,
        )),
        None => Ok(()),
    }
}

/// `HitPosts::write` of query `query`'s report, the reads' texts through `reads`.
fn write_reads(
    report: &[u8],
    hits: &HitPosts,
    query: usize,
    writer: &mut dyn Write,
    sink: &mut dyn Write,
    mut reads: Option<&mut ReadCursor<'_>>,
) -> io::Result<()> {
    hits.write(report, writer, sink, |subject, writer, sink| {
        match reads.as_deref_mut() {
            Some(cursor) => cursor.subject_read(query, subject, writer, sink),
            None => Ok(()),
        }
    })
}

/// The error of a `-html` pairwise report that posts on the diagnostics without the run's warnings
/// (`report_posts` gives them to every such report): LOSAT's internal error, not a silent loss.
pub const HTML_POSTS_WITHOUT_WARNINGS: &str =
    "LOSAT internal error: a -html report posts on the diagnostics, but its format writes no warnings";

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:756-759
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1155-1160
/// Whether the writer of a report of format `format` posts on the run's diagnostics while it
/// writes the report, so that it gets the run's warnings (`QueryWarnings`, whose sink is the
/// diagnostics) whatever its place among the run's formats (the first format gets them anyway, for
/// the posts before each query): BLAST XML and the XML2 and JSON documents (the 803.7 posts and
/// exceptions of their streams, `write_hit_posts`, and the taxonomy database's messages at their
/// `Hit`s), the BLAST archive (NCBI's posts, which it saves as messages) and the Seq-annots of
/// outfmt 8, 9 and 12 (the posts of the bytes that NCBI's stream replaces,
/// `serial::seqalign::write_query_to`), and the pairwise and query-anchored reports (outfmt 0 to 4,
/// `ReportFormat::writes_html_links`) of a run with `-html` (`html`: the posts of their reads of
/// NCBI's formatter registries, `ReportPlace::Html`; line RG T4b fix 1: outfmt 1 to 4 at a later
/// place among a library call's formats posted without them). A writer that posts without
/// warnings is LOSAT's internal error (`write_hit_posts`, `x_pairwise_par::at_place`), not a
/// silent loss.
pub fn report_posts(format: Option<crate::blastinput::app::ReportFormat>, html: bool) -> bool {
    use crate::blastinput::app::ReportFormat;
    matches!(
        format,
        Some(
            ReportFormat::Xml
                | ReportFormat::Xml2
                | ReportFormat::Xml2S
                | ReportFormat::Json
                | ReportFormat::JsonS
                | ReportFormat::ArchiveFormat
                | ReportFormat::AsnText
                | ReportFormat::AsnBinary
                | ReportFormat::JsonSeqalign
        )
    ) || (html && format.is_some_and(ReportFormat::writes_html_links))
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

/// Where the posts of a `-html` report's reads of NCBI's formatter registries go: `posts`, the
/// sink that a caller sets apart for them (`ReportOutputs::html_posts` of a library call's later
/// search, the diagnostics of a later format of BLASTX), else `diagnostics`, the run's. The one
/// owner of that rule, for the reports' writers (`QueryWarnings::html_posts_sink`) and for the
/// prologs' references that post before a search or with its first batch (blastx's
/// `report::write_prolog`, tblastn's and tblastx's prolog writers; line RG T4a2 fix round 3,
/// review e N1).
pub(crate) fn html_posts_sink<'s, W: Write + ?Sized>(
    posts: Option<&'s mut W>,
    diagnostics: &'s mut W,
) -> &'s mut W {
    match posts {
        Some(posts) => posts,
        None => diagnostics,
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:4086
// NCBI reference (598d8ae6): c++/include/corelib/ncbistre.hpp:543-544
/// Writes `text`, what a read of NCBI's formatter registries posts in a `-html` report, to
/// `posts` (`html_posts_sink`); before it, the report written so far (`report`) is flushed where
/// the post goes to the standard error. NCBI's posts go to `cerr` (its default handler), which
/// the C++ library ties to `cout`, so the report reaches standard output first; the posts of a
/// run whose messages follow `-logfile` (tblastx, rpsblast, rpstblastn:
/// `ncbi_environment::run_messages_follow_log`) go to that file and flush nothing, which LOSAT
/// knows by the standard error that the file replaced (`ncbi_environment::terminal_stderr`).
/// NCBI BLAST+ 2.17.0 (line RG T4a2 review e S1, probe `T/t4/review_t4a2_e/out/txlog.rep.*`,
/// `tblastx -html -logfile` with both streams in one pipe): `GETURL_DEBUG`'s `REGISTRY:` line
/// after the version line's flush (166 bytes), where a flush at each post put it after the
/// first description table's head (1,037). A library call's posts go to its diagnostics, never
/// a process stream: its report is flushed before each, as the command line's without
/// `-logfile`.
pub(crate) fn write_html_post<R: Write + ?Sized>(
    report: &mut R,
    posts: &mut dyn Write,
    text: &[u8],
) -> io::Result<()> {
    if text.is_empty() {
        return Ok(());
    }
    if crate::blastinput::ncbi_environment::terminal_stderr().is_none() {
        report.flush()?;
    }
    posts.write_all(text)?;
    posts.flush()
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

impl<'a> QueryWarnings<'a> {
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
    pub fn reads_and_sink(&mut self) -> (Option<&mut ReadCursor<'a>>, &mut dyn Write) {
        (self.reads.as_mut(), &mut *self.sink)
    }

    /// The run ends at the error of a report's writer (NCBI's exception, its stream's failure):
    /// NCBI loads nothing after it, so no text of a read that is not written yet is written, also
    /// not when the warnings end (`ReadCursor::end_run`).
    pub fn end_run(&mut self) {
        if let Some(cursor) = self.reads.as_mut() {
            cursor.end_run();
        }
    }

    /// Writes `text`, what a read of NCBI's formatter registries posts in a `-html` report
    /// (`ReportPlace::Html`), to the sink of the report's posts (`html_posts_sink`), after the
    /// report written so far (`report`) where the post flushes it (`write_html_post`).
    pub fn html_post(&mut self, text: &[u8], report: &mut dyn Write) -> io::Result<()> {
        write_html_post(report, self.html_posts_sink(), text)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3410-3412
    /// Writes `line`, what NCBI's formatter writes to the process's standard error itself in a
    /// `-html` report (`InitConfig`'s `REGISTRY:` line, not a post): with C's `fprintf` to
    /// `stderr`, which does not flush `cout`, so the report written so far is not flushed (in a
    /// pipe of both streams the line comes before the report's bytes since its last flush; line
    /// RG T4a2 fix round 2, review c C4: NCBI BLAST+ 2.17.0 writes it after the prolog's last
    /// `NcbiEndl`, probes `T/t4/fix2/out/c4.*`). It goes to the terminal where the run's messages
    /// follow `-logfile` (tblastx: `ncbi_environment::terminal_stderr`, the standard error before
    /// the command line made it the file; NCBI's posts go to the file, its own line to the
    /// terminal), else to the sink of the report's posts (`html_posts_sink`). A library call's
    /// messages have no terminal of their own: that sink.
    pub fn html_stderr_line(&mut self, line: &[u8]) -> io::Result<()> {
        match crate::blastinput::ncbi_environment::terminal_stderr() {
            Some(mut terminal) => terminal.write_all(line),
            None => {
                let sink = self.html_posts_sink();
                sink.write_all(line)?;
                sink.flush()
            }
        }
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
        sink: &'a mut dyn Write,
    ) -> Self {
        Self {
            before,
            archive,
            sink,
            reads: None,
            later: None,
            html_posts: None,
            report: 0,
        }
    }

    /// These warnings for the report of the search's `report`th format (from 0), whose `-html`
    /// id links call `InitConfig` as a run of that format alone calls it
    /// (`formatter_registry::InitConfig::call_posts`): the first call of each report writes the
    /// call's posts and its `GETURL_DEBUG` line (NCBI's `m_Reg` is a static of the process, which
    /// a run of one format has; line RG T4a2 fix round 2, review c C3 = review d D5, a library
    /// call of formats 0 and 0).
    pub fn for_report(mut self, report: usize) -> Self {
        self.report = report;
        self
    }

    /// These warnings with `posts` (`Some`) as the sink of the posts of a `-html` report's reads
    /// of NCBI's formatter registries (`html_post`, `html_stderr_line`, the prolog's references:
    /// `pairwise::html_with_diagnostics`), where `sink` takes texts that the call writes elsewhere:
    /// the repeats of a later search of a library call (`ReportOutputs::html_posts`) and the
    /// standard error of a later format of BLASTX (`algorithm::blastx::local`), so that each
    /// report posts what a run of its format alone posts, also in a call of several formats
    /// (line RG T4a2 fix round 2, review c C3 = review d D5). `None` leaves them on `sink`.
    pub fn with_html_posts(mut self, posts: Option<&'a mut dyn Write>) -> Self {
        if posts.is_some() {
            self.html_posts = posts;
        }
        self
    }

    /// The report's place among its search's formats (`for_report`).
    pub(crate) fn report(&self) -> usize {
        self.report
    }

    /// The sink of the `-html` report's own posts (`with_html_posts`, else `sink`:
    /// `query_warnings::html_posts_sink`).
    pub(crate) fn html_posts_sink(&mut self) -> &mut dyn Write {
        html_posts_sink(self.html_posts.as_deref_mut(), &mut *self.sink)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
    /// These warnings for a format after the call's first that writes the queries' result
    /// warnings where the first, a structured one, does not: before each query's report the
    /// later lines of `posts[query]` (`Posts::later`) go to `receiver`, the caller's for this
    /// format (`ReportOutputs::later_warnings`), as the command line's run of the format posts
    /// them when the query's report starts. A report that stops has written those of the queries
    /// it reached. Nothing changes without a receiver.
    pub fn with_later(mut self, posts: &'a [Posts], receiver: Option<&'a mut Vec<u8>>) -> Self {
        self.later = receiver.map(|receiver| (posts, receiver));
        self
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:304-308 and 562
    /// Writes the warnings that come after the report of the last of `queries` queries and before
    /// the epilog (psiblast's warnings of the last query's PSSM reads, which it makes after the
    /// query's report): `before[queries]` where the caller sets it, as `before_query`. A report
    /// writer calls it after its last query, before an epilog it writes (review 2 of line DBID
    /// pssmread, H4).
    pub fn after_last_query(&mut self, queries: usize, report: &mut dyn Write) -> io::Result<()> {
        self.before_query(queries, report)
    }

    /// Writes the warnings of query `index`, after the report written so far (`report`)
    /// is flushed; a failed flush returns before them. Then the query's later lines go to the
    /// format's receiver (`with_later`), without a flush: they are no part of the diagnostics.
    pub fn before_query(&mut self, index: usize, report: &mut dyn Write) -> io::Result<()> {
        let placed = before_report(
            self.reads.as_mut(),
            self.before.get(index),
            index,
            report,
            &mut *self.sink,
        )?;
        match placed {
            Some(posts) if posts.has_text() => {
                report.flush()?;
                posts.write(&mut *self.sink, self.archive)?;
                self.sink.flush()
            }
            // A search strategy request's mark alone shows nothing (`Posts::request`).
            Some(posts) => posts.write(&mut *self.sink, self.archive),
            None => Ok(()),
        }?;
        // The query's result warnings of a format after a structured first one (`with_later`).
        if let Some((posts, receiver)) = self.later.as_mut() {
            if let Some(posts) = posts.get(index) {
                posts.write_later(receiver);
            }
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288 (cited at `query_posts`)
/// The writers' step before the report of query `query` (`QueryWarnings::before_query`, and the
/// model of the writers, `LoadTexts::left_before`, which runs it before anything is written): the
/// texts of the reads that the reports of the queries before it did not place
/// (`ReadCursor::subject_read`), which come after those reports (to `sink`, after `report` is
/// flushed), then the posts written before its report (`posts`), whose loads' texts (a batch's
/// prefetch, `Posts::loads`) are made here, where NCBI's prefetch loads the subjects: after the
/// reports of the batches before (`Posts::placed`, in the cursor's ledger). The caller writes the
/// posts returned.
fn before_report(
    cursor: Option<&mut ReadCursor<'_>>,
    posts: Option<&Posts>,
    query: usize,
    report: &mut dyn Write,
    sink: &mut dyn Write,
) -> io::Result<Option<Posts>> {
    let ledger = match cursor {
        Some(cursor) => {
            cursor.before_query(query, report, sink)?;
            cursor.ledger.as_ref()
        }
        None => None,
    };
    Ok(posts.map(|posts| posts.placed_in(ledger)))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:139-158
/// What NCBI writes to the standard error of a query batch's reports before it stops at a row
/// of its last query's report (a `report::tabular::FailedLookupRefusal`), the one owner of that
/// rule for the programs' refusals: for each of the batch's queries up to the row's (`posts`,
/// from the batch's first query; `reads` the same way, shorter where the rest read nothing) the
/// posts before its report (`QueryWarnings::before_query`), the warning of its report's formatter
/// (`formatter`, `FailedLookupRefusal::formatter_warning`), then the texts of the subjects that
/// its rows read one by one (`QueryWarnings::subject_read`, where the batch's prefetch is off);
/// of the row's query only those that its rows up to the row read (`rows`: their subjects, as
/// `reads` names them), which NCBI reads before the row fails.
pub fn posts_before_stop(
    posts: &[Posts],
    reads: &[QueryReads],
    rows: &[u32],
    formatter: Option<&str>,
) -> Posts {
    let mut before = Posts::new();
    let last = posts.len().saturating_sub(1);
    for (query, query_posts) in posts.iter().enumerate() {
        before.append(query_posts);
        if let Some(line) = formatter {
            before.lines(line.as_bytes());
        }
        let query_reads = reads.get(query).map_or(&[][..], Vec::as_slice);
        before.loads(
            query_reads
                .iter()
                .filter(|read| query < last || rows.contains(&read.subject))
                .cloned()
                .collect(),
        );
    }
    before
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:259-262
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288
/// The posts written before each query's report (`QueryWarnings::of_posts`) from the lines of
/// each query (`lines`, each line a warning), NCBI's posts of the bytes replaced in the search
/// strategy requests of the batches (`requests`: each with the batch's first query, numbered as
/// `lines`, in order; `StrategyExport::take_posts`) and the texts of the batches' prefetches
/// (`prefetches`, the same way; `common::database::batch_prefetch`, whose loads' texts are made
/// where they are written, `Posts::loads`). NCBI posts a request's when it writes the request,
/// after it read the batch and before the batch's search, and the prefetch's after the search:
/// in the batch's first query, after the reader's messages of the batch (the first
/// `reader(query)` bytes of its lines), the request's, the prefetch's, then the rest (the query's
/// own warnings). NCBI's posts of its first read of the reported sequences' headers
/// (`header_read`: with the first query of the batch whose Seq-aligns read them,
/// `common::database::Utf8Read`) come after the search, before the prefetch: at the batch's first
/// query after its request's posts and before its prefetch's.
pub fn query_posts(
    lines: Vec<Vec<u8>>,
    requests: Vec<(usize, Posts)>,
    prefetches: Vec<(usize, Posts)>,
    reader: impl Fn(usize) -> usize,
    header_read: Option<(usize, Posts)>,
) -> Vec<Posts> {
    // A batch's first read of a UTF8String header after its request, and its prefetch after both
    // (a stable sort keeps the order of each list).
    let mut requests: Vec<(usize, Posts)> = requests
        .into_iter()
        .chain(header_read)
        .chain(prefetches)
        .collect();
    requests.sort_by_key(|&(first, _)| first);
    let mut requests = requests.into_iter().peekable();
    let mut posts: Vec<Posts> = lines
        .into_iter()
        .enumerate()
        .map(|(query, mut lines)| {
            let Some((_, request)) = requests.next_if(|&(first, _)| first <= query) else {
                return Posts::of_lines(lines);
            };
            let rest = lines.split_off(reader(query).min(lines.len()));
            let mut posts = Posts::of_lines(lines);
            posts.append(&request);
            while let Some((_, request)) = requests.next_if(|&(first, _)| first <= query) {
                posts.append(&request);
            }
            posts.lines(&rest);
            posts
        })
        .collect();
    // A request whose batch starts after the last query's lines (none is made) keeps its posts
    // after them.
    for (_, request) in requests {
        match posts.last_mut() {
            Some(last) => last.append(&request),
            None => posts.push(request),
        }
    }
    posts
}

/// The warning for a query whose ungapped Karlin-Altschul parameters cannot be computed.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2783-2790
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:37-40
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:534-543
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:283-291
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
pub const IDENTITY_CBS_RESET_MESSAGE: &str = "Composition-based statistics cannot be used with the IDENTITY matrix, resetting the composition-based statistics option to 0";

/// NCBI's `kBlastErrMsg_CantCalculateUngappedKAParams` (see `invalid_query_warning`).
pub const INVALID_QUERY_MESSAGE: &str = "Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options";

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:207-214
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1043-1054
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_types.hpp:295-304
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:620-622
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
    warning_line(program, &query_id, messages)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1449-1452
/// The warning line of a query whose id is `query_id`, posted by NCBI's main thread
/// (`ncbi_exception::post_line`: its header alone after the stream failed, line RG D3-1
/// `om/mute/` e01-e04, e12, e16).
fn warning_line(program: &str, query_id: &[u8], messages: &[String]) -> Vec<u8> {
    let mut text = query_id.to_vec();
    text.extend_from_slice(b": ");
    for message in sorted_query_messages(messages) {
        text.extend_from_slice(message.as_bytes());
        text.push(b' ');
    }
    crate::ncbi_exception::post_line(
        crate::ncbi_exception::Severity::Warning,
        None,
        Some(program),
        &text,
    )
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
pub fn query_warning_with_seq_id(
    program: &str,
    seq_id: &[u8],
    title: &[u8],
    messages: &[String],
) -> Vec<u8> {
    if messages.is_empty() {
        return Vec::new();
    }
    warning_line(program, &warning_query_id(seq_id, title), messages)
}

/// The query ID that NCBI's warnings of a query show: the local ID (`Query_N`), then a
/// space and the title when the query has one, cut to its first 25 bytes and `.. ` when
/// it is longer than 35 bytes. The bytes are the title's, raw; the cut can fall inside a
/// UTF-8 sequence.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:533-543
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:920-932
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

/// The error of a protein query with residues that are no NCBIstdaa code (`positions`: their
/// 0-based places in the query), or `None` without any. NCBI's `GetSequenceProtein` throws it;
/// `SetupQueries_OMF` catches it as a warning of the query and invalidates the query's context,
/// and when no context is valid the set-up throws the messages of all queries
/// (`TSearchMessages::ToString`: each one with its severity and a space after it), which the
/// application reports as `BLAST engine error: <that text>` (exit status 3).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:905-917
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-651 and 650-651
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1020-1022 (`ToString`) and
/// c++/include/algo/blast/api/blast_types.hpp:180-184 (`GetMessage` with true)
pub fn invalid_residues_error(positions: &[usize]) -> Option<String> {
    let list = residue_positions(positions)?;
    Some(setup_error(&format!(
        "Invalid residues found at positions {list}"
    )))
}

/// The set-up's error of a query whose descriptors hold an empty reference: NCBI's
/// `CBlastQuerySourceBioseqSet::GetTitle` reads the variant of each, and the null-pointer exception
/// (`ncbi_exception::NULL_POINTER`) is the query's warning, as `invalid_residues_error`'s.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/bioseq_extract_data_priv.cpp:344-372
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:532-543
pub fn null_descriptor_error() -> String {
    setup_error(crate::ncbi_exception::NULL_POINTER)
}

/// The text of the set-up's error of one query whose one message is `message`: its severity
/// (`Warning`) and a space after it (`TSearchMessages::ToString`, cited at
/// `invalid_residues_error`).
fn setup_error(message: &str) -> String {
    format!("Warning: {message} ")
}

/// The list of residue positions of NCBI's `GetSequenceProtein` messages: the first, then up to
/// 19 more after `, `, then `,... (only first 20 shown)` when there are more than 20
/// (`kMaxResiduesToWarnAbout`, blast_setup_cxx.cpp:905-917 above and 920-932, cited at
/// `replaced_o_message`).
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
pub fn few_matches_warning(program: &str) -> Vec<u8> {
    format!("Warning: [{program}] Examining 5 or more matches is recommended\n").into_bytes()
}

#[cfg(test)]
#[path = "query_warnings_tests.rs"]
mod tests;
