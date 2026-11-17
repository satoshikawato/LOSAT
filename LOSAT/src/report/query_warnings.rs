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
    /// The posts to write before the report of each query, by query index (`of_posts`).
    before: &'a [Posts],
    /// The BLAST archive's handle, which saves each whole post as one message (`Posts::write`).
    archive: Option<&'a ArchiveHandle>,
    /// Standard error (the CLI) or the caller's diagnostics.
    pub sink: &'a mut dyn Write,
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1440-1466 (quoted at
// `TaxdbRun::load`): a load asks the taxonomy database for a taxid that its volume's cache does
// not hold, and the run's first ask opens it (the static `CTaxDBFileInfo` of `GetTaxNames`).
/// What the load of a read's subject looks up, as far as it decides the texts: NCBI's lookups write
/// the texts of opening the taxonomy database at the run's first ask, and the offset error at a
/// volume's first lookup of the taxid whose record begins past the end of `taxdb.btd`
/// (`TaxDb::offset_error_taxid`; a volume whose cache would ask again is refused,
/// `TaxdbRun::cache_collision`), nothing else. The model of the reads gives each text to the load
/// that makes the ask, assuming every load it models happens; a load that NCBI does not make (a
/// BLAST XML `Iteration` or `Hit` that the stream leaves out, `ReadCursor::cut`) leaves its texts
/// to NCBI's next load that makes the same ask (`LoadTexts`).
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
}

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
    /// opening and its offset error to the next load that makes those asks (`placed`).
    pub fn left_out(&self, text: &[u8]) {
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

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288 (quoted at
    // `query_posts`: a batch's prefetch, then the result sets of its queries)
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:372-375 (quoted at
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
    ) -> std::sync::Arc<LoadTexts> {
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
                cursor
                    .cut(before, &mut report, &mut sink)
                    .expect("a model writes to no stream");
            }
        }
        cursor
            .before_query(query, &mut report, &mut sink)
            .expect("a model writes to no stream");
        left
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

/// A [`SubjectReadHook`] that writes nothing (the tests of the writers).
#[cfg(test)]
pub fn no_subject_reads<W>(_subject: u32, _writer: &mut W) -> io::Result<()> {
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
    // ```c++
    //     const CSeq_id& kSeqId = kAlignSet.Get().front()->GetSeq_id(1);
    //
    //     try {
    //         const CBioseq_Handle& kSubjBioseqHandle = scope->GetBioseqHandle(kSeqId);
    // ```
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
                (Some(load), Some(ledger)) => load.in_ledger(ledger).left_out(&read.text),
                (Some(load), None) => load.left_out(&read.text),
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

    /// Whether nothing of these posts goes to the diagnostics (no text, no document post).
    fn posts_nothing(&self) -> bool {
        self.texts.is_empty() && self.document.is_empty()
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

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:775-776
// ```c++
// 		s_FillBlastOutput(bxmlout, data);
// 		s_WriteXML2Object(bxmlout, &out_stream);
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
// ```c++
//     catch (const CException& e) {                                           \
//         LOG_POST(Error << "Error: " << e.what());                           \
//         exit_code = BLAST_UNKNOWN_ERROR;                                    \
//     }                                                                       \
// ```
/// The one writer of a structured report's posts (`write_hit_posts`; BLASTX, whose reports keep
/// their own cursor, calls it with that): `report` to `writer` with the posts of its `Hit`s and of
/// its document (`HitPosts::write`) to `sink`, the reads' texts through `reads`. A report that made
/// no read after its last `Hit` (`DocumentPosts::cut`) moves the texts of the later reads
/// (`ReadCursor::cut`); an exception of NCBI's XML stream that ends the run is returned after the
/// posts (`blast_xml::writer::run_ending_exception`), and no text of a read is written after it.
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
            cursor.cut(query, writer, sink)?;
        }
    }
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

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:756-758
// ```c++
// 	catch(CException &e){
// 	    ERR_POST(Error << "Failed s_BlastXMLAddIteration " << query_label << e.what() );
// 	    return;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1155-1160
// ```c++
// void CBlastAppDiagHandler::Post(const SDiagMessage & mess)
// {
// 	if(m_handler != NULL) {
// 		m_handler->Post(mess);
// 	}
// 	if(m_save) {
// ```
/// Whether the writer of a report of format `format` posts on the run's diagnostics while it
/// writes the report, so that it gets the run's warnings (`QueryWarnings`, whose sink is the
/// diagnostics) whatever its place among the run's formats (the first format gets them anyway, for
/// the posts before each query): BLAST XML and the XML2 and JSON documents (the 803.7 posts and
/// exceptions of their streams, `write_hit_posts`, and the taxonomy database's messages at their
/// `Hit`s), the BLAST archive (NCBI's posts, which it saves as messages) and the Seq-annots of
/// outfmt 8, 9 and 12 (the posts of the bytes that NCBI's stream replaces,
/// `serial::seqalign::write_query_to`). A writer that posts without warnings is LOSAT's internal
/// error (`write_hit_posts`), not a silent loss.
pub fn report_posts(format: Option<crate::blastinput::app::ReportFormat>) -> bool {
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
    )
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
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288 (quoted at `query_posts`)
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

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:259-262
// ```c++
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//
//             SaveSearchStrategy(args, m_CmdLineArgs, queries, m_OptsHndl);
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288
// ```c++
//                 BlastFormatter_PreFetchSequenceData(*results, scope,
//                 		                            fmt_args->GetFormattedOutputChoice());
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
// ```
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
    warning_line(program, &query_id, messages)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1449-1452
// ```c++
//     if (!results.GetWarningStrings().empty()) {
//         ...
//         ERR_POST(Warning << results.GetWarningStrings());
// ```
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
    warning_line(program, &warning_query_id(seq_id, title), messages)
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
                SubjectRead::of_text(5, b"A\n"),
                SubjectRead::of_text(7, b"B\n"),
                SubjectRead::of_text(8, b"C\n"),
            ],
            vec![SubjectRead::of_text(3, b"D\n")],
            vec![SubjectRead::of_text(4, b"E\n")],
        ];
        let before: Vec<Posts> = [b"W0\n", b"W1\n", b"W2\n"]
            .map(|lines| Posts::of_lines(lines.to_vec()))
            .into();
        let log = Log::default();
        let (mut report, mut sink) = (log.clone(), log.clone());
        let mut warnings =
            QueryWarnings::of_posts(&before, None, &mut sink).with_subject_reads(&reads);
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
        let reads: Vec<QueryReads> = vec![vec![SubjectRead::of_text(5, b"A\n")]];
        let log = Log::default();
        let (mut report, mut sink) = (log.clone(), log.clone());
        let mut warnings = QueryWarnings::of_posts(&[], None, &mut sink).with_subject_reads(&reads);
        let mut posts = HitPosts::new(warnings.reads_with_texts(0));
        posts.hit(0, 5, |texts| texts.write_all(b"P5\n")).unwrap();
        posts.hit(4, 6, |texts| texts.write_all(b"P6\n")).unwrap();
        write_hit_posts(b"<H5><H6>", &posts, 0, &mut report, Some(&mut warnings)).unwrap();
        drop(warnings);
        let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
        assert_eq!(text, "A\nP5\n<H5>P6\n<H6>");
        // Without warnings, texts that would be lost are LOSAT's internal error; marks alone
        // write the report.
        let mut alone = Vec::new();
        let lost = write_hit_posts(b"<H5><H6>", &posts, 0, &mut alone, None).unwrap_err();
        assert!(
            lost.to_string().starts_with("LOSAT internal error"),
            "{lost}"
        );
        let mut marks = HitPosts::new(true);
        marks.hit(0, 5, |_| Ok(())).unwrap();
        let mut alone = Vec::new();
        write_hit_posts(b"<H5><H6>", &marks, 0, &mut alone, None).unwrap();
        assert_eq!(alone, b"<H5><H6>");
    }

    /// A read with what its load looks up (`ReadLoad`), the run's texts `texts`.
    fn loaded(
        texts: &std::sync::Arc<LoadTexts>,
        subject: u32,
        text: &[u8],
        (volume, offset, opening): (u32, bool, usize),
    ) -> SubjectRead {
        SubjectRead {
            subject,
            text: text.to_vec(),
            load: Some(ReadLoad {
                texts: std::sync::Arc::clone(texts),
                volume,
                looks_up: true,
                offset,
                lead: 0,
                opening,
                scope_utf8: false,
            }),
        }
    }

    /// A BLAST XML report that made no `Hit` after its last mark (`DocumentPosts::cut`: its
    /// `Iteration` or a `Hit` threw) loads none of its later subjects: the texts of their loads go
    /// to NCBI's next load that makes the same ask, of whichever subject (NCBI 2.17.0, review 1 of
    /// K4h M2, `q3w.faa` with `PRE_FETCH_SEQS_LIMIT=1`: the offset error at query 3's `Hit`, not
    /// after query 2's post; agent3 probes m2.*): the texts of opening at the next load that looks
    /// a taxid up, a volume's offset error at its next load that looks that taxid up. A text that
    /// no later load asks for is never written, and a run that NCBI's XML stream ends writes no
    /// text after its posts.
    #[test]
    fn a_cut_report_leaves_its_texts_to_the_next_load_that_asks() {
        let texts = std::sync::Arc::new(LoadTexts::default());
        // Query 0 reads 5 (the run's first ask: the texts of opening, then volume 0's offset
        // error) and 6 (volume 1's offset error); query 1 reads 9 (volume 1, the offset taxid)
        // and 6 again; query 2 reads 7 (volume 0, another taxid).
        let reads: Vec<QueryReads> = vec![
            vec![
                loaded(&texts, 5, b"O\nE0\n", (0, true, 2)),
                loaded(&texts, 6, b"E1\n", (1, true, 0)),
            ],
            vec![
                loaded(&texts, 9, b"", (1, true, 0)),
                loaded(&texts, 6, b"", (1, true, 0)),
            ],
            vec![loaded(&texts, 7, b"", (0, false, 0))],
        ];
        // Query 0 reads 5, then a `Hit` throws: 6's offset error goes to 9's load.
        let log = Log::default();
        let (mut report, mut sink) = (log.clone(), log.clone());
        let mut cursor = ReadCursor::new(&reads).unwrap();
        cursor.subject_read(0, 5, &mut report, &mut sink).unwrap();
        cursor.cut(0, &mut report, &mut sink).unwrap();
        cursor.before_query(1, &mut report, &mut sink).unwrap();
        assert!(cursor.has_texts(1));
        report.write_all(b"<Q1>").unwrap();
        cursor.subject_read(1, 9, &mut report, &mut sink).unwrap();
        cursor.subject_read(1, 6, &mut report, &mut sink).unwrap();
        report.write_all(b"<H6>").unwrap();
        let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
        assert_eq!(text, "O\nE0\n<Q1>E1\n<H6>");
        // Query 0's `Iteration` throws: its texts of opening and volume 1's offset error go to 9's
        // load, volume 0's offset error waits for a load that looks its taxid up (7 does not).
        let log = Log::default();
        let (mut report, mut sink) = (log.clone(), log.clone());
        let mut cursor = ReadCursor::new(&reads).unwrap();
        cursor.cut(0, &mut report, &mut sink).unwrap();
        report.write_all(b"<Q1>").unwrap();
        cursor.subject_read(1, 9, &mut report, &mut sink).unwrap();
        cursor.before_query(2, &mut report, &mut sink).unwrap();
        cursor.subject_read(2, 7, &mut report, &mut sink).unwrap();
        assert!(cursor.unwritten().is_empty());
        let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
        assert_eq!(text, "<Q1>O\nE1\n");
        // A batch modelled after that report: its first read of volume 0 that looks the taxid up
        // makes volume 0's ask, whose text then counts as its own (cut again, it moves again).
        let mut load = reads[0][0].load.clone().unwrap();
        load.opening = 0;
        assert_eq!(load.place(b"X\n".to_vec()), b"X\nE0\n");
        assert_eq!(load.opening, 0);
        load.left_out(b"X\nE0\n");
        let next = loaded(&texts, 8, b"", (0, true, 0));
        assert_eq!(next.placed(), b"X\nE0\n");
        texts.clear();
        // The end of the run: nothing after it.
        let log = Log::default();
        let (mut report, mut sink) = (log.clone(), log.clone());
        let mut cursor = ReadCursor::new(&reads).unwrap();
        cursor.cut(0, &mut report, &mut sink).unwrap();
        cursor.end_run();
        cursor.subject_read(1, 9, &mut report, &mut sink).unwrap();
        cursor.before_query(2, &mut report, &mut sink).unwrap();
        assert!(cursor.unwritten().is_empty() && !cursor.has_texts(1));
        assert!(log.0.lock().unwrap().is_empty());
    }

    /// NCBI's posts of the bytes replaced in a batch's search strategy request follow the reader's
    /// messages of the batch in its first query and come before the rest of its lines (NCBI
    /// 2.17.0, probe RG K/k4b/probe2 `blastp_b3_f6`: the O warning of query 3 after the posts of
    /// its batch); a batch without posts keeps its lines, and the posts are whole posts.
    /// A batch whose prefetch is on, modelled before the report of a batch before it that the
    /// stream cut (a program that models every batch before it writes): its loads' texts are made
    /// where the posts are written (`Posts::loads`), after that report, so the texts that the cut
    /// left go to the prefetch's loads that make their asks, in the order of the loads, after the
    /// batch's reader messages and request posts (NCBI 2.17.0, agent4 m2 probe `*.qm.*.mix.f5`:
    /// the offset error of query 2's prefetch after query 1's `Failed s_BlastXMLAddIteration`).
    #[test]
    fn a_prefetch_takes_the_texts_that_a_cut_report_before_it_left() {
        let texts = std::sync::Arc::new(LoadTexts::default());
        // Query 0's batch reads 5 (the texts of opening, then volume 0's offset error) one by one;
        // query 1's batch prefetches 6 (volume 1, its offset error) and 7 (volume 0, the offset
        // taxid), whose model gave 7 no text: 5 asked before it.
        let reads: Vec<QueryReads> = vec![vec![loaded(&texts, 5, b"O\nE0\n", (0, true, 2))]];
        let mut prefetch = Posts::new();
        prefetch.loads(vec![
            loaded(&texts, 6, b"E1\n", (1, true, 0)),
            loaded(&texts, 7, b"", (0, true, 0)),
        ]);
        let mut request = Posts::new();
        request.post(crate::report::diag_log::ERROR, b"P1\n".to_vec());
        let posts = query_posts(
            vec![Vec::new(), b"R1\nW1\n".to_vec()],
            vec![(1, request)],
            vec![(1, prefetch)],
            |query| [0, 3][query],
            None,
        );
        let log = Log::default();
        let mut report = log.clone();
        let mut sink = log.clone();
        let mut warnings =
            QueryWarnings::of_posts(&posts, None, &mut sink).with_subject_reads(&reads);
        warnings.before_query(0, &mut report).unwrap();
        // Query 0's report is cut before its read of 5.
        let (cursor, sink) = warnings.reads_and_sink();
        cursor.unwrap().cut(0, &mut report, sink).unwrap();
        warnings.before_query(1, &mut report).unwrap();
        drop(warnings);
        let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
        // The opening at 6, the first load that looks a taxid up; volume 0's offset error at 7.
        assert_eq!(text, "R1\nP1\nO\nE1\nE0\nW1\n");
        // Without a cut, the model's texts.
        let texts = std::sync::Arc::new(LoadTexts::default());
        let mut prefetch = Posts::new();
        prefetch.loads(vec![
            loaded(&texts, 6, b"E1\n", (1, true, 0)),
            loaded(&texts, 7, b"", (0, true, 0)),
        ]);
        assert_eq!(prefetch.text(), b"E1\n");
        // A run whose lookups write no text keeps the prefetch's lines as they are.
        let mut plain = Posts::new();
        plain.loads(vec![SubjectRead::of_text(6, b"T\n")]);
        assert_eq!(plain, Posts::of_lines(b"T\n".to_vec()));
    }

    /// The ledger that `LoadTexts::left_before` models for a stopping batch's prefetch is the one
    /// that the writers leave where they have written the reports and posts before it (review 5
    /// of line DBID loadord, M1): a cut report's texts, less what the posts and reads after it take
    /// (the opening at the first load that looks a taxid up, a volume's offset error at its next
    /// load of that taxid); a load of the stop placed from it writes what one placed after the
    /// writers writes.
    #[test]
    fn the_texts_left_before_a_batch_are_those_the_writers_leave() {
        // Query 0's report (cut) reads 5 (the texts of opening, volume 0's offset error) and 6
        // (volume 1's); query 1's batch prefetches 7 (volume 1, the offset taxid: it takes the
        // opening and volume 1's error); query 2's report reads 8 (volume 0, another taxid).
        let shapes = |texts: &std::sync::Arc<LoadTexts>| {
            let reads: Vec<QueryReads> = vec![
                vec![
                    loaded(texts, 5, b"O\nE0\n", (0, true, 2)),
                    loaded(texts, 6, b"E1\n", (1, true, 0)),
                ],
                Vec::new(),
                vec![loaded(texts, 8, b"", (0, false, 0))],
            ];
            let mut prefetch = Posts::new();
            prefetch.loads(vec![loaded(texts, 7, b"", (1, true, 0))]);
            (reads, prefetch)
        };
        let pending = |texts: &LoadTexts| {
            let pending = texts.pending.lock().unwrap();
            (pending.opening.clone(), pending.offsets.clone())
        };
        // The writers: query 0 cut, then query 1's posts, then query 2's report.
        let texts = std::sync::Arc::new(LoadTexts::default());
        let (reads, prefetch) = shapes(&texts);
        let posts = query_posts(
            vec![Vec::new(), Vec::new(), Vec::new()],
            Vec::new(),
            vec![(1, prefetch.clone())],
            |_| 0,
            None,
        );
        let log = Log::default();
        let (mut report, mut sink) = (log.clone(), log.clone());
        let mut warnings =
            QueryWarnings::of_posts(&posts, None, &mut sink).with_subject_reads(&reads);
        warnings.before_query(0, &mut report).unwrap();
        let (cursor, sink) = warnings.reads_and_sink();
        cursor.unwrap().cut(0, &mut report, sink).unwrap();
        warnings.before_query(1, &mut report).unwrap();
        warnings.subject_read(2, 8, &mut report).unwrap();
        warnings.before_query(3, &mut report).unwrap();
        drop(warnings);
        assert_eq!(log.0.lock().unwrap().as_slice(), b"O\nE1\n");
        let written = pending(&texts);
        assert_eq!(written, (Vec::new(), vec![(0, b"E0\n".to_vec())]));
        // The model of the same, before any of it is written, in a run's own ledger.
        let model = std::sync::Arc::new(LoadTexts::default());
        let (reads, prefetch) = shapes(&model);
        let prefetched = |query: usize| (query == 1).then_some(&prefetch);
        let left = LoadTexts::left_before(&reads, &|query| query == 0, &prefetched, 3);
        assert_eq!(pending(&left), written);
        assert!(pending(&model) == (Vec::new(), Vec::new()) && !left.is_empty());
        // The stop's load of volume 0 that looks the taxid up writes volume 0's error.
        let mut stop = loaded(&model, 9, b"", (0, true, 0));
        stop.place_left(&left);
        assert_eq!(stop.text, b"E0\n");
        assert!(left.is_empty());
        // Up to query 1 only (a stop in query 1's batch): its prefetch is that batch's own.
        let left = LoadTexts::left_before(&reads, &|query| query == 0, &prefetched, 1);
        assert_eq!(
            pending(&left),
            (
                b"O\n".to_vec(),
                vec![(0, b"E0\n".to_vec()), (1, b"E1\n".to_vec())]
            )
        );
        // No cut: nothing left.
        assert!(LoadTexts::left_before(&reads, &|_| false, &prefetched, 3).is_empty());
    }

    #[test]
    fn request_posts_follow_the_reader_messages_of_their_batch() {
        use crate::report::diag_log::ERROR;
        let lines = vec![
            b"R0\nW0\n".to_vec(),
            b"W1\n".to_vec(),
            b"R2\nT2\nW2\n".to_vec(),
        ];
        let post = |text: &[u8]| {
            let mut posts = Posts::new();
            posts.post(ERROR, text);
            posts
        };
        let reader = [3, 0, 3];
        let posts = query_posts(
            lines,
            vec![(0, post(b"P0\n")), (2, post(b"P2a\n")), (2, post(b"P2b\n"))],
            Vec::new(),
            |query| reader[query],
            None,
        );
        let texts: Vec<Vec<u8>> = posts.iter().map(Posts::text).collect();
        assert_eq!(
            texts,
            [
                b"R0\nP0\nW0\n".to_vec(),
                b"W1\n".to_vec(),
                b"R2\nP2a\nP2b\nT2\nW2\n".to_vec()
            ]
        );
        let mut whole = Posts::of_lines(b"R2\n".to_vec());
        whole.append(&post(b"P2a\n"));
        whole.append(&post(b"P2b\n"));
        whole.lines(b"T2\nW2\n");
        assert_eq!(posts[2], whole);
        // A request after the last query's lines keeps its posts after them.
        let posts = query_posts(
            vec![b"W0\n".to_vec()],
            vec![(1, post(b"P1\n"))],
            Vec::new(),
            |_| 0,
            None,
        );
        assert_eq!(posts[0].text(), b"W0\nP1\n");
        // NCBI's first read of a UTF8String header (RG-T5d fix 1) after the request's posts of
        // its batch, before the rest; a batch without a request has it after its reader messages.
        let header = || Some((2, Posts::of_lines(b"U2\n".to_vec())));
        let lines = || vec![b"R0\nW0\n".to_vec(), b"W1\n".to_vec(), b"R2\nT2\n".to_vec()];
        let with = query_posts(
            lines(),
            vec![(0, post(b"P0\n")), (2, post(b"P2\n"))],
            vec![(2, Posts::of_lines(b"F2\n".to_vec()))],
            |q| reader[q],
            header(),
        );
        assert_eq!(with[2].text(), b"R2\nP2\nU2\nF2\nT2\n");
        let without = query_posts(lines(), Vec::new(), Vec::new(), |q| reader[q], header());
        assert_eq!(without[2].text(), b"R2\nU2\nT2\n");
        assert_eq!(without[0].text(), b"R0\nW0\n");
    }

    #[test]
    fn a_query_without_texts_reads_nothing() {
        let reads: Vec<QueryReads> = vec![Vec::new(), Vec::new()];
        let mut sink = Vec::new();
        let warnings = QueryWarnings::of_posts(&[], None, &mut sink).with_subject_reads(&reads);
        assert!(!warnings.reads_with_texts(0));
        assert!(ReadCursor::new(&reads).is_none());
    }
}
