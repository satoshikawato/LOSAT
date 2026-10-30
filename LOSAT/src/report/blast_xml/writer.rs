//! The type writers of the `CBlastOutput` model (objects/blastxml/blastxml.asn) over an
//! `Encoder`, and BLAST's incremental document (`BlastXML_FormatReport`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:725-733
//! ```c++
//!     if( incremental_output && incremental_struct->m_IterationNum == 0) {
//!         bool add_dtd_reference = true, add_xml_version = true;
//!         s_SerializeAndSplitBy( bxmlout, "</BlastOutput_param>",
//!                 serial_xml_start, serial_xml_end,
//!                 add_dtd_reference, add_xml_version );
//!         // incremental_output
//!         *out_stream << serial_xml_start << "\n<BlastOutput_iterations>" ;
//!         incremental_struct->m_SerialXmlEnd = "\n</BlastOutput_iterations>" + serial_xml_end;
//!     }
//! ```
//! So a document is: once, the `BlastOutput` object (filled from the first formatted query)
//! cut after `</BlastOutput_param>` and `\n<BlastOutput_iterations>`; per query the
//! `Iteration` object cut after `</Iteration_query-len>`, `\n<Iteration_hits>\n`, each `Hit`
//! written alone, `</Iteration_hits>` and the `Iteration` tail (`Statistics`, `message`,
//! `\n</Iteration>`); at the end `\n</BlastOutput_iterations>\n</BlastOutput>\n` and the
//! epilog's end of line (`CBlastFormat::PrintEpilog`, blast_format.cpp:2242-2243).

use anyhow::Result;

use crate::report::serial::encoder::{Encoder, SerialError};
use crate::report::serial::xml::XmlEncoder;

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
        // `if ((val = data->GetMatchReward()) != 0) params.SetSc_match(val);`
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
    enc.end_sequence();
    enc.end_member();
}

/// The `BlastOutput` object through `</BlastOutput_param>`, filled from the first formatted
/// query (`BlastXML_FormatReport`, blastxml_format.cpp:679-720), then `\n<BlastOutput_iterations>`.
fn write_header(out: &mut Vec<u8>, run: &XmlRun<'_>, first: &XmlQuery<'_>) -> Option<SerialError> {
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
    enc.error()
}

/// `Statistics` (`s_BlastXMLGetStatistics`, blastxml_format.cpp:628-650).
fn write_statistics<E: Encoder>(enc: &mut E, run: &XmlRun<'_>, stats: &XmlStats) {
    enc.member("stat", 5);
    enc.begin_sequence("Statistics");
    int_member(enc, "db-num", 0, run.db_num);
    enc.member("db-len", 1);
    enc.big_int(run.db_len);
    enc.end_member();
    int_member(enc, "hsp-len", 2, stats.hsp_len);
    // `stat->SetEff_space((double)data->GetEffectiveSearchSpace(index))`
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
    pub(crate) midline: Vec<u8>,
    pub(crate) columns: Vec<usize>,
    pub(crate) segments: Vec<(usize, usize)>,
    hit: Vec<u8>,
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
}

/// The hit-level calls of a structured writer behind `HitSink`.
pub trait HitWrite {
    fn hit(
        &mut self,
        hit: &XmlHit<'_>,
        hsps: &mut dyn FnMut(&mut HspSink<'_>) -> Result<()>,
    ) -> Result<()>;
}

/// The HSP-level call of a structured writer behind `HspSink`.
pub trait HspWrite {
    fn hsp(&mut self, hsp: &XmlHsp<'_>) -> Result<()>;
}

/// The state of a document: the number of queries formatted so far (NCBI's
/// `SBlastXMLIncremental::m_IterationNum`), which numbers the iterations and tells whether the
/// header is still to be written.
pub struct BlastXmlWriter {
    iterations: u32,
    scratch: Scratch,
}

impl BlastXmlWriter {
    /// A writer that has formatted `queries_before` queries (0 for a document's start; the
    /// count of the earlier calls when a program formats a document in several calls).
    pub fn new(queries_before: usize) -> Self {
        Self {
            iterations: queries_before as u32,
            scratch: Scratch::default(),
        }
    }
}

impl StructuredXml for BlastXmlWriter {
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
        if self.iterations == 0 {
            if let Some(error) = write_header(out, run, query) {
                return Err(error.into());
            }
        }
        self.iterations += 1;
        let mut enc = XmlEncoder::object(out);
        enc.begin_sequence("Iteration");
        int_member(&mut enc, "iter-num", 0, self.iterations as i32);
        string_member(&mut enc, "query-ID", 1, query.id);
        string_member(&mut enc, "query-def", 2, query_def(query.title));
        int_member(&mut enc, "query-len", 3, query.len);
        // `Iteration_hits`: the Hit objects are written one by one between the head and the tail
        // (s_BlastXMLAddIteration, blastxml_format.cpp:576-589).
        enc.raw(b"\n<Iteration_hits>\n");
        let mut writer = Xml5Hits {
            run,
            masks: query.masks,
            iteration: &mut enc,
            scratch: &mut self.scratch,
            hits: 0,
        };
        hits(&mut HitSink::new(&mut writer))?;
        enc.raw(b"</Iteration_hits>");
        write_statistics(&mut enc, run, &query.stats);
        if let Some(message) = query.message {
            string_member(&mut enc, "message", 6, message);
        }
        enc.end_sequence();
        match enc.error() {
            Some(error) => Err(error.into()),
            None => Ok(()),
        }
    }

    /// Appends the end of the document if any query was formatted:
    /// `\n</BlastOutput_iterations>\n</BlastOutput>\n` and the epilog's end of line. With no
    /// query formatted the end of the document is empty and only the epilog's end of line
    /// remains.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2242-2243
    /// ```c++
    ///     if (m_FormatType == CFormattingArgs::eXml) {
    ///         m_Outfile << m_BlastXMLIncremental->m_SerialXmlEnd << endl;
    /// ```
    /// (`m_SerialXmlEnd` is set by the first formatted query, blastxml_format.cpp:755-763.)
    fn finish(&mut self, out: &mut Vec<u8>) -> Result<()> {
        if self.iterations > 0 {
            out.extend_from_slice(b"\n</BlastOutput_iterations>\n</BlastOutput>\n\n");
        } else {
            out.push(b'\n');
        }
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
}

/// The hit level of outfmt 5: one `Hit` object per subject, written alone and joined into the
/// iteration.
struct Xml5Hits<'a, 'o> {
    run: &'a XmlRun<'a>,
    masks: &'a [QueryMask],
    iteration: &'a mut XmlEncoder<'o>,
    scratch: &'a mut Scratch,
    hits: i32,
}

impl HitWrite for Xml5Hits<'_, '_> {
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
            int_member(&mut enc, "len", 4, hit.len);
            enc.member("hsps", 5);
            enc.begin_list(false);
            let mut writer = Xml5Hsps {
                run: self.run,
                masks: self.masks,
                hit: &mut enc,
                scratch,
                hsps: 0,
            };
            hsps(&mut HspSink::new(&mut writer))?;
            enc.end_list();
            enc.end_member();
            enc.end_sequence();
            enc.end_document();
            match enc.error() {
                Some(error) => Err(error.into()),
                None => Ok(()),
            }
        })();
        if result.is_ok() {
            self.iteration.raw(&hit_buf);
        }
        scratch.hit = hit_buf;
        result
    }
}

/// The HSP level of outfmt 5.
struct Xml5Hsps<'a, 'o> {
    run: &'a XmlRun<'a>,
    masks: &'a [QueryMask],
    hit: &'a mut XmlEncoder<'o>,
    scratch: &'a mut Scratch,
    hsps: i32,
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
/// midline (`scratch.midline`) and the three counts.
pub(crate) struct HspCounts {
    pub identity: i32,
    pub positive: i32,
    pub gaps: i32,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:238-241
    // ```c++
    //         	if(!kTranslated)
    //         	{
    //         		num_gaps = kAlign.GetTotalGapCount();
    // ```
    // The gaps of the Seq-align's segments (also for a translated search, lines 259-260): the
    // columns of the edit script's gaps, not the `-` letters of the rows, which a stored gap
    // character of a protein subject also is. Without an edit script, the `-` letters.
    let gaps = match hsp.ops {
        Some(ops) => ops
            .iter()
            .filter(|op| !op.is_sub())
            .map(|op| op.num() as usize)
            .sum::<usize>() as i32,
        None => hsp
            .qrow
            .iter()
            .chain(hsp.srow)
            .filter(|&&letter| letter == b'-')
            .count() as i32,
    };
    HspCounts {
        identity,
        positive,
        gaps,
    }
}

impl HspWrite for Xml5Hsps<'_, '_> {
    fn hsp(&mut self, hsp: &XmlHsp<'_>) -> Result<()> {
        self.hsps += 1;
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
        string_member(enc, "hseq", 18, hsp.srow);
        string_member(enc, "midline", 19, &scratch.midline);
        enc.end_sequence();
        Ok(())
    }
}

#[cfg(test)]
mod segment_tests {
    use super::*;

    #[test]
    fn an_empty_edit_script_is_one_segment_like_no_script() {
        let mut segments = Vec::new();
        for ops in [None, Some(&[][..])] {
            fill_segments(&mut segments, ops, 7);
            assert_eq!(segments, [(0, 6)]);
        }
        fill_segments(&mut segments, Some(&[][..]), 0);
        assert!(segments.is_empty());
        let ops = [GapEditOp::Sub(3), GapEditOp::Ins(2), GapEditOp::Sub(4)];
        fill_segments(&mut segments, Some(&ops), 9);
        assert_eq!(segments, [(0, 2), (3, 4), (5, 8)]);
    }
}
