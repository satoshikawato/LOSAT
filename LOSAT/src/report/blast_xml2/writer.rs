//! The type writers of the `CBlastOutput2` model (objects/blastxml2/blastxml2.asn) over an
//! `Encoder`, one report per query, and `Xml2Writer`, which frames the reports (`framing`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:678-686
//! A report is the `BlastOutput2` object of one query: the program, version and reference, the
//! target (the `-db` name or the `Subject_N` list), the parameters and the results (one
//! `Search`; a `-subject` search's is the one element of `bl2seq`).

use anyhow::Result;

use crate::report::blast_xml::hsp::{round6, MaskLetter};
use crate::report::blast_xml::records::{HitTaxonomy, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats};
use crate::report::blast_xml::writer::{
    count_hsp, int_member, real_member, replaced_posts, string_member, DocumentPost, HitSink,
    HitWrite, HspSink, HspWrite, Scratch, StructuredXml,
};
use crate::report::serial::encoder::{Encoder, SerialError};
use crate::report::serial::json::JsonEncoder;
use crate::report::serial::xml::XmlEncoder;

use super::framing::{MadeReport, NoSideFiles, RunSideFiles, SideFiles, Xml2Document, Xml2Kind};
use super::records::{Xml2Run, Xml2Target};

/// `Parameters` (`s_SetBlastXMLParameters`, blastxml2_format.cpp:514-555), members in the order
/// of blastxml2.asn: the matrix, `expect`, the nucleotide scores, the gap costs of a gapped
/// search, the filter, the composition-based statistics and the genetic codes.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:531-535
fn write_parameters<E: Encoder>(enc: &mut E, run: &XmlRun<'_>, params: &Xml2Run) {
    enc.member("params", 4);
    enc.begin_sequence("Parameters");
    if let Some(matrix) = run.matrix {
        string_member(enc, "matrix", 0, matrix.as_bytes());
    }
    real_member(enc, "expect", 1, run.expect);
    if let Some((reward, penalty)) = run.reward_penalty {
        if reward != 0 {
            int_member(enc, "sc-match", 3, reward);
        }
        if penalty != 0 {
            int_member(enc, "sc-mismatch", 4, penalty);
        }
    }
    if params.gapped {
        int_member(enc, "gap-open", 5, run.gap_open);
        int_member(enc, "gap-extend", 6, run.gap_extend);
    }
    if !run.filter.is_empty() {
        string_member(enc, "filter", 7, run.filter);
    }
    if params.cbs != 0 {
        int_member(enc, "cbs", 10, params.cbs);
    }
    if let Some(code) = params.query_gencode {
        int_member(enc, "query-gencode", 11, code);
    }
    if let Some(code) = params.db_gencode {
        int_member(enc, "db-gencode", 12, code);
    }
    enc.end_sequence();
    enc.end_member();
}

/// `Statistics` (`s_SetBlastXMLStatistics`, blastxml2_format.cpp:558-577): the database's
/// sequence and letter counts only for `-db`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:565-568
fn write_statistics<E: Encoder>(enc: &mut E, run: &XmlRun<'_>, params: &Xml2Run, stats: &XmlStats) {
    enc.member("stat", 5);
    enc.begin_sequence("Statistics");
    if params.target == Xml2Target::Database {
        enc.member("db-num", 0);
        enc.big_int(run.db_num);
        enc.end_member();
        enc.member("db-len", 1);
        enc.big_int(run.db_len);
        enc.end_member();
    }
    int_member(enc, "hsp-len", 2, stats.hsp_len);
    // `eff-space` is a BigInt here (a REAL in outfmt 5).
    enc.member("eff-space", 3);
    enc.big_int(stats.eff_space);
    enc.end_member();
    real_member(enc, "kappa", 4, stats.kappa);
    real_member(enc, "lambda", 5, stats.lambda);
    real_member(enc, "entropy", 6, stats.entropy);
    enc.end_sequence();
    enc.end_member();
}

/// `NCBI_BlastOutput2::kNoHitsFound` joined to the message of a query of a batch that NCBI did
/// not search (blastxml2_format.cpp:614-622).
const UNSEARCHED_SUFFIX: &[u8] = b" \nNo hits found";

/// One `Search`: the query, the hits of `hits`, the statistics and the message
/// (`s_SetBlastXMLSearch`, blastxml2_format.cpp:579-633).
fn write_search<E, F>(
    enc: &mut E,
    run: &XmlRun<'_>,
    params: &Xml2Run,
    query: &XmlQuery<'_>,
    scratch: &mut Scratch,
    start: usize,
    hits: F,
) -> Result<()>
where
    E: Encoder,
    F: FnOnce(&mut HitSink<'_>) -> Result<()>,
{
    enc.begin_sequence("Search");
    string_member(enc, "query-id", 0, query.id);
    // NCBI sets the query title only when it is not empty.
    if !query.title.is_empty() {
        string_member(enc, "query-title", 1, query.title);
    }
    int_member(enc, "query-len", 2, query.len);
    // The masks of a query that has any, but not those of a minus frame. A query of a batch that
    // was not searched lists none: NCBI makes its result without the search, so the masks the
    // search would have reported (probe: `BATCH_SIZE=10 -lcase_masking`, an N-only query with a
    // lower-case run) are not there.
    if query.searched && !query.masks.is_empty() {
        enc.member("query-masking", 3);
        enc.begin_list(false);
        for mask in query.masks.iter().filter(|mask| mask.frame >= 0) {
            enc.element();
            enc.begin_sequence("Range");
            int_member(enc, "from", 0, mask.from);
            int_member(enc, "to", 1, mask.to);
            enc.end_sequence();
        }
        enc.end_list();
        enc.end_member();
    }
    // A query of a batch that was not searched has no `hits` member.
    if query.searched {
        enc.member("hits", 4);
        enc.begin_list(false);
        {
            let mut writer = Xml2Hits {
                enc: &mut *enc,
                run,
                params,
                masks: query.masks,
                scratch: &mut *scratch,
                hits: 0,
                start,
            };
            hits(&mut HitSink::new(&mut writer))?;
        }
        enc.end_list();
        enc.end_member();
    }
    write_statistics(enc, run, params, &query.stats);
    if query.searched {
        if let Some(message) = query.message {
            string_member(enc, "message", 6, message);
        }
    } else {
        let mut message = query.message.unwrap_or_default().to_vec();
        message.extend_from_slice(UNSEARCHED_SUFFIX);
        string_member(enc, "message", 6, &message);
    }
    enc.end_sequence();
    Ok(())
}

/// The `BlastOutput2` of one query.
fn write_report<E, F>(
    enc: &mut E,
    run: &XmlRun<'_>,
    params: &Xml2Run,
    query: &XmlQuery<'_>,
    scratch: &mut Scratch,
    start: usize,
    hits: F,
) -> Result<()>
where
    E: Encoder,
    F: FnOnce(&mut HitSink<'_>) -> Result<()>,
{
    enc.begin_document("BlastOutput2");
    enc.begin_sequence("BlastOutput2");
    enc.member("report", 0);
    enc.begin_sequence("Report");
    string_member(enc, "program", 0, run.program.as_bytes());
    string_member(enc, "version", 1, run.version.as_bytes());
    string_member(enc, "reference", 2, run.reference.as_bytes());
    enc.member("search-target", 3);
    match params.target {
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xml2format.cpp:157-169
        // `i->name` is each name of the `-db` list as the reports open it alone
        // (`seqdb::alias::report_names`): quotes kept, one space between names.
        Xml2Target::Database => {
            enc.variant("Target", "db", 0);
            // The names' bytes reach the stream as they are; the stream writes or throws on them.
            let names = crate::seqdb::alias::report_names(run.db_name).join(" ");
            enc.visible_str(names.as_bytes());
            enc.end_variant();
        }
        Xml2Target::Subjects(count) => match &params.subject_names {
            Some(names) => {
                enc.variant_string_list("Target", "subjects", 1, names.iter().cloned());
            }
            None => {
                enc.variant_string_list(
                    "Target",
                    "subjects",
                    1,
                    (1..=count).map(|n| format!("Subject_{n}")),
                );
            }
        },
    }
    enc.end_member();
    write_parameters(enc, run, params);
    enc.member("results", 5);
    match params.target {
        Xml2Target::Database => {
            enc.variant("Results", "search", 1);
            write_search(enc, run, params, query, scratch, start, hits)?;
            enc.end_variant();
        }
        Xml2Target::Subjects(_) => {
            enc.variant("Results", "bl2seq", 2);
            enc.begin_list(false);
            enc.element();
            write_search(enc, run, params, query, scratch, start, hits)?;
            enc.end_list();
            enc.end_variant();
        }
    }
    enc.end_member();
    enc.end_sequence();
    enc.end_member();
    enc.end_sequence();
    enc.end_document();
    match enc.error() {
        Some(error) => Err(error.into()),
        None => Ok(()),
    }
}

/// The hit level of an XML2 `Search` (`s_SeqAlignToXMLHit`, blastxml2_format.cpp:347-432).
struct Xml2Hits<'a, E> {
    enc: &'a mut E,
    run: &'a XmlRun<'a>,
    params: &'a Xml2Run,
    masks: &'a [crate::report::blast_xml::QueryMask],
    scratch: &'a mut Scratch,
    hits: i32,
    /// Where the report that this call writes starts in the output (`Xml2Document::report`):
    /// the query's own or, for a `-subject` run, the one held before it (no database, so no
    /// posts at its `Hit`s).
    start: usize,
}

impl<E: Encoder> HitWrite for Xml2Hits<'_, E> {
    /// The start of the report written now: NCBI makes the whole report's data, loading the
    /// subject of each `Hit`, before it writes any of it (`Xml2Document::report`).
    fn position(&self) -> usize {
        self.start
    }

    fn hit(
        &mut self,
        hit: &XmlHit<'_>,
        hsps: &mut dyn FnMut(&mut HspSink<'_>) -> Result<()>,
    ) -> Result<()> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:412-417 (cited
        // at `report::db_subject::ParsedSubject::lookup_throws`)
        // The `catch` of the `Hit` whose subject's lookup throws loads the subject again outside
        // any `catch` of the hit: the exception leaves the query's report (`CaughtException`),
        // before anything of the `Hit` is made.
        if let Some(exception) = hit.lookup_throws {
            return Err(crate::report::blast_xml::CaughtException(exception.clone()).into());
        }
        self.hits += 1;
        let database = self.params.target == Xml2Target::Database;
        let enc = &mut *self.enc;
        enc.element();
        enc.begin_sequence("Hit");
        int_member(enc, "num", 0, self.hits);
        // One `HitDescr` per Blast-def-line: LOSAT's databases have one per sequence, a
        // `-subject` record one.
        enc.member("description", 1);
        enc.begin_list(false);
        // The `HitDescr` of each definition line of a sequence with parsed Seq-ids
        // (`report::db_subject::ParsedSubject::hit_descriptions`).
        for description in hit.descriptions {
            enc.element();
            enc.begin_sequence("HitDescr");
            string_member(enc, "id", 0, &description.id);
            if !description.accession.is_empty() {
                string_member(enc, "accession", 1, &description.accession);
            }
            if let Some(title) = &description.title {
                string_member(enc, "title", 2, title);
            }
            if description.taxid != 0 {
                int_member(enc, "taxid", 3, description.taxid);
                if let HitTaxonomy::Named(sciname) = &description.taxonomy {
                    string_member(enc, "sciname", 4, sciname);
                }
            }
            enc.end_sequence();
        }
        if hit.descriptions.is_empty() {
            enc.element();
            enc.begin_sequence("HitDescr");
            string_member(enc, "id", 0, hit.id);
            // A taxid the taxonomy database does not know (`GetTaxInfo` throws) leaves the id alone
            // (blastxml2_format.cpp:422-428).
            if !matches!(hit.taxonomy, HitTaxonomy::NotFound(_)) {
                // `-subject` ids are local: no accession (blastxml2_format.cpp:398-402).
                if database {
                    string_member(enc, "accession", 1, hit.accession);
                }
                string_member(enc, "title", 2, hit.title);
                if database && hit.taxid != 0 {
                    int_member(enc, "taxid", 3, hit.taxid);
                    if let HitTaxonomy::Named(sciname) = hit.taxonomy {
                        string_member(enc, "sciname", 4, sciname);
                    }
                }
            }
            enc.end_sequence();
        }
        enc.end_list();
        enc.end_member();
        int_member(enc, "len", 2, hit.len);
        enc.member("hsps", 3);
        enc.begin_list(false);
        {
            let mut writer = Xml2Hsps {
                enc: &mut *enc,
                run: self.run,
                masks: self.masks,
                scratch: &mut *self.scratch,
                hsps: 0,
            };
            hsps(&mut HspSink::new(&mut writer))?;
        }
        enc.end_list();
        enc.end_member();
        enc.end_sequence();
        Ok(())
    }
}

/// The HSP level of an XML2 `Search` (`s_SeqAlignSetToXMLHsps`, blastxml2_format.cpp:95-339).
struct Xml2Hsps<'a, E> {
    enc: &'a mut E,
    run: &'a XmlRun<'a>,
    masks: &'a [crate::report::blast_xml::QueryMask],
    scratch: &'a mut Scratch,
    hsps: i32,
}

impl<E: Encoder> HspWrite for Xml2Hsps<'_, E> {
    fn hsp(&mut self, hsp: &XmlHsp<'_>) -> Result<()> {
        self.hsps += 1;
        let program = self.run.program;
        let blastn = program == "blastn";
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:159-170
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:247-259
        // A frame is written for a nucleotide sequence of a Std-seg (translated) Seq-align: the
        // query of blastx, tblastx and rpstblastn (whose profiles are proteins), the subject of
        // tblastn and tblastx.
        let query_translated = matches!(program, "blastx" | "tblastx" | "rpstblastn");
        let subject_translated = matches!(program, "tblastn" | "tblastx");
        // The query's regions are applied by the frame of a translated query, else frame 0
        // (blastxml2_format.cpp:211, :253): a nucleotide query's regions of frame 0 are
        // applied too, which outfmt 5 does not do.
        let mask_frame = if query_translated { hsp.query_frame } else { 0 };
        let counts = count_hsp(
            hsp,
            self.run.positives,
            self.masks,
            mask_frame,
            MaskLetter::Lower,
            &mut *self.scratch,
        );
        let enc = &mut *self.enc;
        enc.element();
        enc.begin_sequence("Hsp");
        int_member(enc, "num", 0, self.hsps);
        real_member(enc, "bit-score", 1, round6(hsp.bit_score));
        real_member(enc, "score", 2, f64::from(hsp.score));
        real_member(enc, "evalue", 3, round6(hsp.evalue));
        int_member(enc, "identity", 4, counts.identity);
        if !blastn {
            int_member(enc, "positive", 5, counts.positive);
        }
        int_member(enc, "query-from", 9, hsp.query_from);
        int_member(enc, "query-to", 10, hsp.query_to);
        if blastn {
            string_member(enc, "query-strand", 11, b"Plus");
        }
        if query_translated {
            int_member(enc, "query-frame", 12, hsp.query_frame);
        }
        int_member(enc, "hit-from", 13, hsp.hit_from);
        int_member(enc, "hit-to", 14, hsp.hit_to);
        if blastn {
            let strand: &[u8] = if hsp.hit_from > hsp.hit_to {
                b"Minus"
            } else {
                b"Plus"
            };
            string_member(enc, "hit-strand", 15, strand);
        }
        if subject_translated {
            int_member(enc, "hit-frame", 16, hsp.hit_frame);
        }
        int_member(enc, "align-len", 17, hsp.qrow.len() as i32);
        int_member(enc, "gaps", 18, counts.gaps);
        string_member(enc, "qseq", 19, &self.scratch.masked);
        string_member(enc, "hseq", 20, counts.hseq(hsp.srow, self.scratch));
        string_member(enc, "midline", 21, &self.scratch.midline);
        enc.end_sequence();
        Ok(())
    }
}

/// Whether the text of the element `tag` of a report is a VisibleString (`string_member`, the
/// target's `db` and `subjects`) rather than a number (`int_member`, `real_member`, `big_int`):
/// NCBI's stream writes a string one character at a time and a number at once, which decides
/// where its buffer reaches the output (`out_buffer::xml_written_at`).
pub(super) fn is_string_element(tag: &[u8]) -> bool {
    matches!(
        tag,
        b"program"
            | b"version"
            | b"reference"
            | b"db"
            | b"subjects"
            | b"matrix"
            | b"filter"
            | b"query-id"
            | b"query-title"
            | b"id"
            | b"accession"
            | b"title"
            | b"sciname"
            | b"query-strand"
            | b"hit-strand"
            | b"qseq"
            | b"hseq"
            | b"midline"
            | b"message"
    )
}

/// What an XML report of outfmt 16 or 14 hands the framing (`MadeReport`) after its writer's
/// `result`: the posts of its replaced control characters and, when NCBI's stream threw, the
/// report cut at the failing byte (`XmlThrow::at`) and the exception's report. Outfmt 16 catches
/// the exception and posts it after the bytes written (the stream's destructor writes its buffer
/// out), and the run goes on; outfmt 14 does not catch it (`MadeReport::ends_run`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:749-766
/// The report of 16 is written to the output stream through the object stream's buffer
/// (`s_WriteXML2ObjectNoHeader`, blastxml2_format.cpp:716-724): a post comes after the bytes that
/// the buffer had written by then (`out_buffer::xml_written_at`). A 14 report goes to its file,
/// so its posts only keep their order.
fn xml_report_made(
    report: &mut Vec<u8>,
    kind: Xml2Kind,
    program: &str,
    result: Result<()>,
    posts: Vec<crate::report::serial::xml::XmlPost>,
    thrown: Option<Box<crate::report::serial::xml::XmlThrow>>,
) -> Result<MadeReport> {
    match (result, thrown) {
        (Err(error), thrown)
            if thrown.is_none()
                || !matches!(
                    error.downcast_ref::<SerialError>(),
                    Some(SerialError::Utf8(_))
                ) =>
        {
            Err(error)
        }
        (_, thrown) => {
            if let Some(thrown) = &thrown {
                report.truncate(thrown.at);
            }
            let written: Vec<u64> = if kind == Xml2Kind::Single16 && !posts.is_empty() {
                let at: Vec<usize> = posts.iter().map(|post| post.at).collect();
                crate::report::serial::out_buffer::xml_written_at(report, &at, is_string_element)
            } else {
                vec![0; posts.len()]
            };
            let mut written = written.into_iter();
            let mut made = MadeReport {
                posts: replaced_posts(posts, program, |_| written.next().unwrap_or(0) as usize),
                ends_run: None,
            };
            match (kind, thrown) {
                (Xml2Kind::Single16, Some(thrown)) => {
                    let exception = thrown.exception();
                    let message = exception.what_after(exception.get_msg());
                    made.posts.push(DocumentPost {
                        at: report.len(),
                        text: crate::ncbi_exception::error_post_line(None, Some(program), &message),
                    });
                }
                (_, thrown) => made.ends_run = thrown,
            }
            Ok(made)
        }
    }
}

/// The writer of the XML2 and JSON documents of outfmt 16, 15, 14 and 13: one complete report per
/// query, framed by `Xml2Document`.
pub struct Xml2Writer<F: SideFiles = NoSideFiles> {
    document: Xml2Document<F>,
    params: Xml2Run,
    scratch: Scratch,
    report: Vec<u8>,
}

impl Xml2Writer<NoSideFiles> {
    /// A writer of the document `kind` for a search of `params.target`; a `-subject` search
    /// holds the last report back one query (`Xml2Document`). The unit tests make writers this
    /// way; a run uses `for_run` or `with_document`.
    #[cfg(test)]
    pub fn new(kind: Xml2Kind, params: Xml2Run) -> Self {
        let delayed = params.target.holds_back();
        Self::with_document(Xml2Document::new(kind, delayed), params)
    }
}

impl Xml2Writer<RunSideFiles> {
    /// A writer of the document `kind` for the run with `-out` `out`: outfmt 13 and 14 write
    /// their side files next to it, or to the caller's receiver `sink`, and the master file to
    /// the output stream (`Xml2Document::for_run`).
    pub fn for_run(
        kind: Xml2Kind,
        params: Xml2Run,
        out: Option<&std::path::Path>,
        sink: Option<&super::SideFileSink>,
    ) -> Self {
        let delayed = params.target.holds_back();
        Self::with_document(Xml2Document::for_run(kind, delayed, out, sink), params)
    }
}

impl<F: SideFiles> Xml2Writer<F> {
    pub fn with_document(document: Xml2Document<F>, params: Xml2Run) -> Self {
        Self {
            document,
            params,
            scratch: Scratch::default(),
            report: Vec::new(),
        }
    }

    /// The framing of the writer's document, for a run whose reports are written by several
    /// writers in turn (`Xml2Document::take_carry` and `resume`).
    pub fn document_mut(&mut self) -> &mut Xml2Document<F> {
        &mut self.document
    }
}

impl<F: SideFiles> StructuredXml for Xml2Writer<F> {
    fn query<H>(
        &mut self,
        run: &XmlRun<'_>,
        query: &XmlQuery<'_>,
        out: &mut Vec<u8>,
        hits: H,
    ) -> Result<()>
    where
        H: FnOnce(&mut HitSink<'_>) -> Result<()>,
    {
        let (params, scratch) = (&self.params, &mut self.scratch);
        let kind = self.document.kind();
        // The report is made after the framing's separator, where the `Hit`s load their subjects
        // (`HitSink::position`: the start of the report written now).
        self.document
            .report(&mut self.report, out, |start, report| match kind {
                Xml2Kind::Single16 | Xml2Kind::Files14 => {
                    let mut enc = if kind == Xml2Kind::Single16 {
                        XmlEncoder::xml2_report(report)
                    } else {
                        XmlEncoder::xml2_document(report)
                    };
                    let result = write_report(&mut enc, run, params, query, scratch, start, hits);
                    let posts = enc.take_posts();
                    let thrown = enc.take_thrown();
                    xml_report_made(report, kind, run.program, result, posts, thrown)
                }
                Xml2Kind::Single15 => {
                    let mut enc = JsonEncoder::array_element(report);
                    write_report(&mut enc, run, params, query, scratch, start, hits)
                        .map(|()| MadeReport::default())
                }
                Xml2Kind::Files13 => {
                    let mut enc = JsonEncoder::new(report);
                    write_report(&mut enc, run, params, query, scratch, start, hits)
                        .map(|()| MadeReport::default())
                }
            })
    }

    fn take_posts(&mut self) -> crate::report::blast_xml::DocumentPosts {
        self.document.take_posts()
    }

    fn start(&mut self, reported: Option<&XmlRun<'_>>, out: &mut Vec<u8>) -> Result<()> {
        if let Some(run) = reported {
            run.check_sums(match self.document.kind() {
                Xml2Kind::Files13 => 13,
                Xml2Kind::Files14 => 14,
                Xml2Kind::Single15 => 15,
                Xml2Kind::Single16 => 16,
            })?;
        }
        self.document.start(out);
        Ok(())
    }

    fn finish(&mut self, out: &mut Vec<u8>) -> Result<()> {
        self.document.epilog(out)
    }

    fn holds_back(&self) -> bool {
        self.document.holds_back()
    }

    fn unwritten(&mut self, out: &mut Vec<u8>) -> Result<()> {
        self.document.unwritten(out)
    }
}
