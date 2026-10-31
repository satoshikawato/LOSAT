//! SAM (outfmt 17) of the gapped Seq-aligns of a nucleotide search: NCBI's `CBlast_SAM_Formatter`
//! (`CSAM_Formatter`, `CSAM_CIGAR_Formatter`, `CCIGAR_Formatter`), which keeps every record until
//! the formatter is destroyed and writes the header last. One [`SamDocument`] holds a run's
//! records as rendered bytes, its `@SQ` lines and the references it has listed; [`SamDocument::finish`]
//! writes the document.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:313-316, :405-442
//! ```c++
//! CSAM_Formatter::~CSAM_Formatter(void)
//! {
//!     Flush();
//! }
//! ...
//! void CSAM_Formatter::Flush(void)
//! {
//!     if ( !m_Out ) return;
//!     // Headers
//!     bool have_data = !m_Header.m_Data.empty() || !m_Body.empty();
//!     if (have_data) {
//!         *m_Out << "@HD\tVN:1.2";
//!         x_PrintSOTag();
//!         x_PrintGOTag();
//!         *m_Out << '\n';
//!     }
//!     ITERATE(CSAM_Headers::TData, it, m_Header.m_Data) {
//!         *m_Out << it->second << '\n';
//!     }
//!     if (have_data  &&  !m_ProgramInfo.m_Id.empty()) {
//!         *m_Out << "@PG\tID:" << m_ProgramInfo.m_Id;
//!         ...
//!     }
//!     // Alignments
//!     ITERATE(TLines, it, m_Body) {
//!         *m_Out << *it << '\n';
//!     }
//!     m_Header.m_Data.clear();
//!     m_Body.clear();
//! }
//! ```
//!
//! The record of an HSP is the one NCBI's `CSAM_CIGAR_Formatter::AddRow` makes from the HSP's
//! Dense-seg (`CCIGAR_Formatter::x_FormatLine` for the CIGAR); the Dense-seg is the one that
//! outfmt 8 and 12 write (`report::serial::seqalign::gapped_segments`).

use std::collections::HashSet;
use std::io::{self, Write};

use anyhow::Result;

use super::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use super::serial::seqalign::{gapped_segments, seqalign_evalue, Segments, Strand, GAP_VALUE};

/// The words of the custom specification of outfmt 17 that NCBI reads.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/sam.cpp:84-104
/// ```c++
/// void CBlast_SAM_Formatter::x_ProcessCustomSpec(const string & custom_spec,
///                                                const CSAM_Formatter::SProgramInfo & info)
/// {
///     vector<string> format_tokens;
///     NStr::Split(custom_spec, " ", format_tokens);
///     CSAM_Formatter::SetProgram(info);
///     m_refRow = 1;
///     ITERATE (vector<string>, iter, format_tokens) {
///          if("SR" == *iter) {
///             m_refRow = 0;
///          }
///          if ("SQ" == *iter) {
///             CSAM_Formatter::SetFlag(CSAM_Formatter::fSAM_SeqData);
///         }
///     }
///
///     if(m_refRow == 1) {
///         CSAM_Formatter::SetGroupOrder(eGO_Reference);
///         CSAM_Formatter::SetSortOrder(eSO_Coordinate);
///     }
/// }
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SamSpec {
    /// `SR`: the subject is the reference and the query the target (`m_refRow` 0).
    pub subject_as_reference: bool,
    /// `SQ`: the target's letters are written in SEQ (`fSAM_SeqData`).
    pub sequence: bool,
}

impl SamSpec {
    /// The words of `spec` (split at every space; exactly `SQ` and `SR` count, anything else
    /// is ignored).
    pub fn parse(spec: &str) -> Self {
        let mut parsed = Self::default();
        for word in spec.split(' ') {
            match word {
                "SR" => parsed.subject_as_reference = true,
                "SQ" => parsed.sequence = true,
                _ => {}
            }
        }
        parsed
    }
}

/// A sequence that a record names: the query, or a subject.
#[derive(Clone, Copy, Debug)]
pub struct SamSequence<'a> {
    /// The Seq-id and the searched length (`QueryRecord`, `SubjectRecord`) that the Dense-seg
    /// is made from.
    pub record: SamId<'a>,
    /// A number that is the same for the same sequence in every call of a run: the references
    /// listed in `@SQ` are told apart by it.
    pub key: u64,
    /// The sequence's index for [`SamLetters`].
    pub index: usize,
    /// The length of the whole Bioseq (`GetBioseqLength`, also where a range was searched).
    pub length: i32,
}

/// The Seq-id, searched length and range start of a query or a subject.
#[derive(Clone, Copy, Debug)]
pub enum SamId<'a> {
    Query(QueryRecord<'a>),
    Subject(SubjectRecord<'a>),
}

impl SamId<'_> {
    fn name(&self, out: &mut Vec<u8>) {
        match self {
            SamId::Query(query) => query.id.sam_name(out),
            SamId::Subject(subject) => subject.id.sam_name(out),
        }
    }
}

/// Which row's sequence the letters are of.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SamRow {
    Query,
    Subject,
}

/// The letters of the Bioseqs for SEQ.
pub trait SamLetters {
    /// Appends the letters `[from, to)` (record coordinates) of query or subject `index` to
    /// `out`, in upper case IUPAC.
    fn letters(
        &mut self,
        row: SamRow,
        index: usize,
        from: usize,
        to: usize,
        out: &mut Vec<u8>,
    ) -> Result<()>;
}

/// No letters (the specification has no `SQ`).
pub struct NoLetters;

impl SamLetters for NoLetters {
    fn letters(&mut self, _: SamRow, _: usize, _: usize, _: usize, _: &mut Vec<u8>) -> Result<()> {
        Ok(())
    }
}

/// Bytes held in blocks of one MiB, so a long document never copies what it holds.
#[derive(Default)]
struct ChunkedBytes {
    blocks: Vec<Vec<u8>>,
}

impl ChunkedBytes {
    const BLOCK: usize = 1 << 20;

    fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }

    fn push(&mut self, bytes: &[u8]) {
        match self.blocks.last_mut() {
            Some(last) if last.capacity() - last.len() >= bytes.len() => {
                last.extend_from_slice(bytes)
            }
            _ => {
                let mut block = Vec::with_capacity(Self::BLOCK.max(bytes.len()));
                block.extend_from_slice(bytes);
                self.blocks.push(block);
            }
        }
    }
}

/// A rendered record of one query, before the sort of the default mode.
struct Rendered {
    /// `min(GetSeqStart(0), GetSeqStop(0))` of the query row, then the e-value
    /// (`SortHspByMasterStartAscending`).
    start: i32,
    evalue: f64,
    /// The bytes of the record in `Scratch::buf`.
    bytes: std::ops::Range<usize>,
    /// The subject of the record (index into the call's subjects).
    subject: usize,
}

/// Buffers reused from HSP to HSP and from query to query.
#[derive(Default)]
struct Scratch {
    segments: Segments,
    rows: Vec<Rendered>,
    buf: Vec<u8>,
    letters: Vec<u8>,
    name: Vec<u8>,
    ops: Vec<(u8, i64)>,
}

/// The document of a run.
pub struct SamDocument {
    spec: SamSpec,
    /// The `@PG` line without its CL tag and newline: `@PG\tID:0\tVN:<version>`.
    program: Vec<u8>,
    /// The raw command line: every argument followed by one space; `None`: no CL tag.
    command_line: Option<Vec<u8>>,
    /// The `@SQ` lines, one per reference in order of first appearance.
    sq: Vec<u8>,
    seen: HashSet<u64>,
    body: ChunkedBytes,
    scratch: Scratch,
}

impl SamDocument {
    /// A document for `spec`, the program `blastn` of version `version` (`2.17.0+`), and the
    /// command line that the `@PG` line shows (`None`: NCBI's formatter without one writes no
    /// CL tag; an empty line is the same).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2460-2466
    /// ```c++
    ///     CSAM_Formatter::SProgramInfo info("0", CBlastVersion().Print(), m_Cmdline);
    ///     info.m_Name = m_Program;
    /// ```
    pub fn new(spec: SamSpec, version: &str, command_line: Option<Vec<u8>>) -> Self {
        Self {
            spec,
            program: format!("@PG\tID:0\tVN:{version}").into_bytes(),
            command_line: command_line.filter(|line| !line.is_empty()),
            sq: Vec::new(),
            seen: HashSet::new(),
            body: ChunkedBytes::default(),
            scratch: Scratch::default(),
        }
    }

    /// Adds the records of one query's final hit list: `hsps`, grouped by subject in the
    /// Seq-align order, `subjects` indexed by `HspRecord::subject`. `run` says whether the
    /// subjects are remapped (`-subject`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/sam.cpp:50-61
    /// ```c++
    /// void CBlast_SAM_Formatter::x_Print(const CSeq_align_set & aln)
    /// {
    ///     if(m_refRow == 1) {
    ///         CSeq_align_set sorted;
    ///         sorted.Set() = aln.Get();
    ///         sorted.Set().sort(align_format::CAlignFormatUtil::SortHspByMasterStartAscending);
    ///         CSAM_Formatter::Print(sorted, m_refRow);
    ///     }
    ///     else {
    ///         CSAM_Formatter::Print(aln, m_refRow);
    ///     }
    /// }
    /// ```
    pub fn query(
        &mut self,
        run: &RunRecord<'_>,
        query: &SamSequence<'_>,
        subjects: &[SamSequence<'_>],
        hsps: &[HspRecord<'_>],
        letters: &mut dyn SamLetters,
    ) -> Result<()> {
        if hsps.is_empty() {
            return Ok(());
        }
        let mut scratch = std::mem::take(&mut self.scratch);
        scratch.rows.clear();
        scratch.buf.clear();
        let result = (|| -> Result<()> {
            for hsp in hsps {
                let subject = &subjects[hsp.subject as usize];
                self.render(run, query, subject, hsp, letters, &mut scratch)?;
                scratch.rows.last_mut().expect("a row was rendered").subject = hsp.subject as usize;
            }
            if !self.spec.subject_as_reference {
                // `list::sort` is stable.
                scratch.rows.sort_by(|a, b| {
                    a.start.cmp(&b.start).then(
                        a.evalue
                            .partial_cmp(&b.evalue)
                            .unwrap_or(std::cmp::Ordering::Equal),
                    )
                });
            }
            for row in &scratch.rows {
                let reference = if self.spec.subject_as_reference {
                    &subjects[row.subject]
                } else {
                    query
                };
                self.list_reference(reference, &mut scratch.name);
                self.body.push(&scratch.buf[row.bytes.clone()]);
            }
            Ok(())
        })();
        // A huge query (a chromosome against a chromosome) renders its records here before they
        // are copied into the body: its buffers do not keep that capacity for the rest of the run.
        const KEEP_BYTES: usize = 4 << 20;
        const KEEP_ROWS: usize = 1 << 16;
        if scratch.buf.capacity() > KEEP_BYTES {
            scratch.buf.clear();
            scratch.buf.shrink_to(KEEP_BYTES / 4);
        }
        if scratch.letters.capacity() > KEEP_BYTES {
            scratch.letters.clear();
            scratch.letters.shrink_to(KEEP_BYTES / 4);
        }
        if scratch.rows.capacity() > KEEP_ROWS {
            scratch.rows.clear();
            scratch.rows.shrink_to(KEEP_ROWS / 4);
        }
        self.scratch = scratch;
        result
    }

    /// Lists `reference` in `@SQ` the first time it is seen.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:181-187, :445-452
    /// ```c++
    ///     CBioseq_Handle refseq = GetScope()->GetBioseqHandle(GetRefId());
    ///     if (m_KnownRefSeqs.find(refseq) == m_KnownRefSeqs.end()) {
    ///         m_Head.AddSequence(CSeq_id_Handle::GetHandle(GetRefId()),
    ///             "@SQ\tSN:" + x_GetRefIdString() +
    ///             "\tLN:" + NStr::UInt8ToString(refseq.GetBioseqLength()));
    ///         m_KnownRefSeqs.insert(refseq);
    ///     }
    /// ```
    fn list_reference(&mut self, reference: &SamSequence<'_>, name: &mut Vec<u8>) {
        if !self.seen.insert(reference.key) {
            return;
        }
        name.clear();
        reference.record.name(name);
        self.sq.extend_from_slice(b"@SQ\tSN:");
        self.sq.extend_from_slice(name);
        self.sq.extend_from_slice(b"\tLN:");
        push_int(&mut self.sq, i64::from(reference.length));
        self.sq.push(b'\n');
    }

    /// Renders the record of one HSP into `scratch` (`CSAM_CIGAR_Formatter::AddRow` after
    /// `CCIGAR_Formatter::x_FormatLine`).
    fn render(
        &self,
        run: &RunRecord<'_>,
        query: &SamSequence<'_>,
        subject: &SamSequence<'_>,
        hsp: &HspRecord<'_>,
        letters: &mut dyn SamLetters,
        scratch: &mut Scratch,
    ) -> Result<()> {
        let (SamId::Query(query_record), SamId::Subject(subject_record)) =
            (&query.record, &subject.record)
        else {
            unreachable!("the first row is the query and the second the subject")
        };
        gapped_segments(
            run,
            query_record,
            subject_record,
            hsp,
            &mut scratch.segments,
        );
        let segments = &scratch.segments;
        let starts = segments.starts();
        let lens = segments.lens();
        let strands = segments.strands();
        let offsets = [query_record.offset, subject_record.offset];
        let start_of = |seg: usize, row: usize| -> i32 {
            let start = starts[2 * seg + row];
            if start >= 0 {
                start + offsets[row]
            } else {
                start
            }
        };
        // The reference row and the target row (`FormatByTargetRow(m_refRow)`).
        let subject_ref = self.spec.subject_as_reference;
        let (tgt_row, ref_row) = if subject_ref { (0, 1) } else { (1, 0) };
        // `StrandSign(row)`: the row's strand in the first segment.
        let reverse = (strands[tgt_row] == Strand::Minus) != (strands[ref_row] == Strand::Minus);
        let ref_minus = strands[ref_row] == Strand::Minus;
        let numseg = lens.len();
        // The CIGAR (`CCIGAR_Formatter::x_FormatLine`, widths 1), merging neighbours of the
        // same type; the ranges of the rows are the unions of their segments.
        let ops = &mut scratch.ops;
        ops.clear();
        let mut ref_range = (i64::MAX, i64::MIN);
        let mut tgt_range = (i64::MAX, i64::MIN);
        let mut aligned: i64 = 0;
        for idx in 0..numseg {
            let seg = if ref_minus { numseg - idx - 1 } else { idx };
            let len = i64::from(lens[seg]);
            let ref_start = start_of(seg, ref_row);
            let tgt_start = start_of(seg, tgt_row);
            let ref_seq = ref_start != GAP_VALUE;
            let tgt_seq = tgt_start != GAP_VALUE;
            let kind = if tgt_seq && !ref_seq {
                b'I'
            } else if !tgt_seq && ref_seq {
                b'D'
            } else if tgt_seq && ref_seq {
                aligned += len;
                b'M'
            } else {
                continue;
            };
            if ref_seq {
                ref_range.0 = ref_range.0.min(i64::from(ref_start));
                ref_range.1 = ref_range.1.max(i64::from(ref_start) + len);
            }
            if tgt_seq {
                tgt_range.0 = tgt_range.0.min(i64::from(tgt_start));
                tgt_range.1 = tgt_range.1.max(i64::from(tgt_start) + len);
            }
            match ops.last_mut() {
                Some((last, count)) if *last == kind => *count += len,
                _ => ops.push((kind, len)),
            }
        }
        let num_dif: i64 = ops
            .iter()
            .filter(|(kind, _)| *kind != b'M')
            .map(|(_, count)| *count)
            .sum();
        // The target's clips (`AddRow`).
        let (target, reference) = if subject_ref {
            (query, subject)
        } else {
            (subject, query)
        };
        let target_length = i64::from(target.length);
        let mut clip_front: Option<i64> = None;
        let mut clip_back: Option<i64> = None;
        if tgt_range.0 > 0 {
            if reverse {
                clip_back = Some(tgt_range.0);
            } else {
                clip_front = Some(tgt_range.0);
            }
        }
        if tgt_range.1 < target_length {
            if reverse {
                clip_front = Some(target_length - tgt_range.1);
            } else {
                clip_back = Some(target_length - tgt_range.1);
            }
        }
        // SEQ (`fSAM_SeqData`).
        let want_letters = self.spec.sequence;
        scratch.letters.clear();
        if want_letters {
            let row = if subject_ref {
                SamRow::Query
            } else {
                SamRow::Subject
            };
            letters.letters(
                row,
                target.index,
                tgt_range.0 as usize,
                tgt_range.1 as usize,
                &mut scratch.letters,
            )?;
            if reverse {
                reverse_complement(&mut scratch.letters);
            }
        }
        // The record.
        let begin = scratch.buf.len();
        let buf = &mut scratch.buf;
        target.record.name(buf);
        buf.push(b'\t');
        buf.extend_from_slice(if reverse { b"16" } else { b"0" });
        buf.push(b'\t');
        reference.record.name(buf);
        buf.push(b'\t');
        push_int(buf, ref_range.0 + 1);
        buf.extend_from_slice(b"\t255\t");
        if let Some(front) = clip_front {
            push_int(buf, front);
            buf.push(b'H');
        }
        for &(kind, count) in ops.iter() {
            push_int(buf, count);
            buf.push(kind);
        }
        if let Some(back) = clip_back {
            push_int(buf, back);
            buf.push(b'H');
        }
        buf.extend_from_slice(b"\t*\t0\t0\t");
        if want_letters {
            buf.extend_from_slice(&scratch.letters);
        } else {
            buf.push(b'*');
        }
        buf.extend_from_slice(b"\t*");
        // The tags: AS, EV, NM, PI, BS, in the order of the code (`AS + EV + NM + PI + BS`).
        buf.extend_from_slice(b"\tAS:i:");
        push_int(buf, i64::from(hsp.score));
        buf.extend_from_slice(b"\tEV:f:");
        let evalue = seqalign_evalue(hsp.evalue);
        push_g(buf, evalue);
        buf.extend_from_slice(b"\tNM:i:");
        push_int(buf, num_dif);
        buf.extend_from_slice(b"\tPI:f:");
        let identity = if i64::from(hsp.num_ident) == aligned {
            100.0
        } else {
            f64::min(99.99, 100.0 * f64::from(hsp.num_ident) / aligned as f64)
        };
        push_fixed2(buf, identity);
        buf.extend_from_slice(b"\tBS:f:");
        push_g(buf, hsp.bit_score);
        buf.push(b'\n');
        let end = buf.len();
        scratch.rows.push(Rendered {
            start: query_row_start(starts, offsets[0]),
            evalue,
            bytes: begin..end,
            subject: 0,
        });
        Ok(())
    }

    /// Writes the document and empties it: nothing when no record was added; else `@HD`, the
    /// `@SQ` lines, `@PG` and the records. A second call writes nothing.
    pub fn finish(&mut self, out: &mut dyn Write) -> io::Result<()> {
        if self.body.is_empty() {
            return Ok(());
        }
        if self.spec.subject_as_reference {
            out.write_all(b"@HD\tVN:1.2\tGO:query\n")?;
        } else {
            out.write_all(b"@HD\tVN:1.2\tSO:coordinate\tGO:reference\n")?;
        }
        out.write_all(&std::mem::take(&mut self.sq))?;
        out.write_all(&self.program)?;
        if let Some(command_line) = self.command_line.as_deref() {
            out.write_all(b"\tCL:")?;
            out.write_all(command_line)?;
        }
        out.write_all(b"\tPN:blastn\n")?;
        for block in std::mem::take(&mut self.body.blocks) {
            out.write_all(&block)?;
        }
        self.seen.clear();
        Ok(())
    }

    /// Whether the document holds a record.
    pub fn has_records(&self) -> bool {
        !self.body.is_empty()
    }
}

/// The lowest coordinate of the query row in record coordinates:
/// `min(GetSeqStart(0), GetSeqStop(0))` of the Seq-align (the minimum over the segments that
/// are not gaps).
fn query_row_start(starts: &[i32], offset: i32) -> i32 {
    starts
        .chunks_exact(2)
        .map(|segment| segment[0])
        .filter(|&start| start >= 0)
        .min()
        .map_or(0, |start| start + offset)
}

/// Appends the decimal digits of `value`.
fn push_int(out: &mut Vec<u8>, value: i64) {
    out.extend_from_slice(value.to_string().as_bytes());
}

/// The reverse complement of upper case IUPAC letters, in place (`eNa_strand_minus` of
/// `CSeqVector` in `eCoding_Iupac`).
pub fn reverse_complement(letters: &mut [u8]) {
    letters.reverse();
    for letter in letters {
        *letter = match *letter {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            other => other,
        };
    }
}

/// `%.2f` (`NStr::DoubleToString(value, 2)`).
fn push_fixed2(out: &mut Vec<u8>, value: f64) {
    out.extend_from_slice(format!("{value:.2}").as_bytes());
}

/// C's `%g` (precision 6), which `NStr::DoubleToString(value)` prints.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2316-2340
/// ```c++
/// void NStr::DoubleToString(string& out_str, double value,
///                           int precision, TNumToStringFlags flags)
/// {
///     ...
///     if (precision >= 0 || ...) {
///         ...
///     } else {
///         const char* format;
///         switch (flags & fDoubleGeneral) {
///             ...
///             case fDoubleGeneral: // default
///             default:
///                 format = "%g";
///                 break;
///         }
///         n = ::snprintf(buffer, kMaxDoubleStringSize, format, value);
/// ```
pub fn push_g(out: &mut Vec<u8>, value: f64) {
    if value.is_nan() {
        out.extend_from_slice(if value.is_sign_negative() {
            b"-nan"
        } else {
            b"nan"
        });
        return;
    }
    if value.is_sign_negative() {
        out.push(b'-');
    }
    let value = value.abs();
    if value.is_infinite() {
        out.extend_from_slice(b"inf");
        return;
    }
    if value == 0.0 {
        out.push(b'0');
        return;
    }
    // The exponent after rounding to 6 significant digits decides the style.
    let scientific = format!("{value:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').expect("an exponent");
    let exponent: i32 = exponent.parse().expect("a decimal exponent");
    if !(-4..6).contains(&exponent) {
        let mantissa = if mantissa.contains('.') {
            mantissa.trim_end_matches('0').trim_end_matches('.')
        } else {
            mantissa
        };
        out.extend_from_slice(mantissa.as_bytes());
        out.push(b'e');
        out.push(if exponent < 0 { b'-' } else { b'+' });
        let magnitude = exponent.unsigned_abs();
        if magnitude < 10 {
            out.push(b'0');
        }
        out.extend_from_slice(magnitude.to_string().as_bytes());
    } else {
        let places = (5 - exponent) as usize;
        let fixed = format!("{value:.places$}");
        let fixed = if fixed.contains('.') {
            fixed.trim_end_matches('0').trim_end_matches('.')
        } else {
            &fixed
        };
        out.extend_from_slice(fixed.as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::GapEditOp;
    use crate::report::seqid::SeqId;

    fn g(value: f64) -> String {
        let mut out = Vec::new();
        push_g(&mut out, value);
        String::from_utf8(out).unwrap()
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2316-2340 (`%g`).
    #[test]
    fn g_prints_six_significant_digits_as_c_does() {
        assert_eq!(g(0.0), "0");
        assert_eq!(g(1e-180), "1e-180");
        assert_eq!(g(2.60916e-178), "2.60916e-178");
        assert_eq!(g(0.8216932), "0.821693");
        assert_eq!(g(608.6682), "608.668");
        assert_eq!(g(55395.12), "55395.1");
        assert_eq!(g(100000.0), "100000");
        assert_eq!(g(999999.5), "1e+06");
        assert_eq!(g(1e6), "1e+06");
        assert_eq!(g(123456789.0), "1.23457e+08");
        assert_eq!(g(0.0001), "0.0001");
        assert_eq!(g(0.00001), "1e-05");
        assert_eq!(g(0.5), "0.5");
        assert_eq!(g(2.0), "2");
        assert_eq!(g(1.5e-5), "1.5e-05");
        assert_eq!(g(-0.25), "-0.25");
        assert_eq!(g(f64::INFINITY), "inf");
    }

    #[test]
    fn the_specification_words_are_exact() {
        let parse = SamSpec::parse;
        assert_eq!(parse(""), SamSpec::default());
        assert_eq!(
            parse("SQ"),
            SamSpec {
                subject_as_reference: false,
                sequence: true
            }
        );
        assert_eq!(
            parse("SR SQ"),
            SamSpec {
                subject_as_reference: true,
                sequence: true
            }
        );
        for plain in ["sq", "foo SQ,SR", "std", "SQSR", "S R"] {
            assert_eq!(parse(plain), SamSpec::default(), "{plain}");
        }
        assert_eq!(parse("foo SQ").sequence, true);
        assert_eq!(parse("SQ  SQ").subject_as_reference, false);
    }

    #[test]
    fn a_reverse_complement_keeps_the_ambiguity_codes() {
        let mut letters = b"ACGTRYKMBVDHNSW".to_vec();
        reverse_complement(&mut letters);
        assert_eq!(letters, b"WSNDHBVKMRYACGT");
    }

    /// A document of one query with the given HSP over a subject of length 100; the query
    /// has 60 letters.
    fn document(
        spec: SamSpec,
        cl: Option<&[u8]>,
        hsp: HspRecord<'_>,
        letters: &mut dyn SamLetters,
    ) -> String {
        let run = RunRecord {
            program: crate::blastinput::app::Program::Blastn,
            db_title: None,
            gapped: true,
        };
        let query = SamSequence {
            record: SamId::Query(QueryRecord {
                id: SeqId::Local("Query_1"),
                searched_length: 60,
                offset: 0,
            }),
            key: 0,
            index: 0,
            length: 60,
        };
        let subject = SamSequence {
            record: SamId::Subject(SubjectRecord {
                id: SeqId::Local("Subject_1"),
                searched_length: 100,
                offset: 0,
                use_this_seqid: &[],
            }),
            key: 0,
            index: 0,
            length: 100,
        };
        let mut document = SamDocument::new(spec, "2.17.0+", cl.map(<[u8]>::to_vec));
        document
            .query(&run, &query, &[subject], &[hsp], letters)
            .unwrap();
        let mut out = Vec::new();
        document.finish(&mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    fn hsp(
        q: (i32, i32),
        s: (i32, i32),
        frame: i8,
        ident: i32,
        ops: &[GapEditOp],
    ) -> HspRecord<'_> {
        HspRecord {
            subject: 0,
            q_offset: q.0,
            q_end: q.1,
            s_offset: s.0,
            s_end: s.1,
            q_frame: frame,
            s_frame: 1,
            score: 40,
            num: 0,
            evalue: 2.5e-7,
            bit_score: 75.8,
            num_ident: ident,
            num_positives: 0,
            comp_adjust: 0,
            ops,
        }
    }

    struct Fixed;

    impl SamLetters for Fixed {
        fn letters(
            &mut self,
            row: SamRow,
            _: usize,
            from: usize,
            to: usize,
            out: &mut Vec<u8>,
        ) -> Result<()> {
            let alphabet: &[u8] = if row == SamRow::Query {
                b"AACCGGTTAC"
            } else {
                b"ACGTRYACGT"
            };
            out.extend((from..to).map(|at| alphabet[at % alphabet.len()]));
            Ok(())
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:179-298 (AddRow) and
    // cigar_formatter.cpp:157-352 (x_FormatLine): the target is the subject, the reference the
    // query; a plus-strand HSP with a gap of 3 in the subject (D) and one of 2 in the query (I).
    #[test]
    fn a_plus_hsp_has_the_subject_as_target_and_clips_it() {
        let ops = [
            GapEditOp::Sub(10),
            GapEditOp::Ins(3),
            GapEditOp::Sub(20),
            GapEditOp::Del(2),
            GapEditOp::Sub(5),
        ];
        // The query consumes 10 + 3 + 20 + 5 = 38 letters, the subject 10 + 20 + 2 + 5 = 37.
        let hsp = hsp((4, 42), (20, 57), 1, 30, &ops);
        let text = document(SamSpec::default(), Some(b"blastn -x "), hsp, &mut NoLetters);
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "@HD\tVN:1.2\tSO:coordinate\tGO:reference");
        assert_eq!(lines[1], "@SQ\tSN:Query_1\tLN:60");
        assert_eq!(lines[2], "@PG\tID:0\tVN:2.17.0+\tCL:blastn -x \tPN:blastn");
        let fields: Vec<&str> = lines[3].split('\t').collect();
        assert_eq!(fields[0], "Subject_1");
        assert_eq!(fields[1], "0");
        assert_eq!(fields[2], "Query_1");
        assert_eq!(fields[3], "5");
        assert_eq!(fields[4], "255");
        // The clips come from the subject range [20, 57) of a length of 100.
        assert_eq!(fields[5], "20H10M3D20M2I5M43H");
        assert_eq!(fields[9], "*");
        assert_eq!(fields[11], "AS:i:40");
        assert_eq!(fields[12], "EV:f:2.5e-07");
        assert_eq!(fields[13], "NM:i:5");
        assert_eq!(fields[14], "PI:f:85.71");
        assert_eq!(fields[15], "BS:f:75.8");
        assert_eq!(lines.len(), 4);
    }

    #[test]
    fn sq_writes_the_target_letters_and_sr_swaps_the_roles() {
        let ops = [GapEditOp::Sub(10), GapEditOp::Ins(3), GapEditOp::Sub(10)];
        let spec = SamSpec::parse("SQ SR");
        let hsp = hsp((4, 27), (20, 40), 1, 20, &ops);
        let text = document(spec, None, hsp, &mut Fixed);
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "@HD\tVN:1.2\tGO:query");
        assert_eq!(lines[1], "@SQ\tSN:Subject_1\tLN:100");
        // No CL tag without a command line.
        assert_eq!(lines[2], "@PG\tID:0\tVN:2.17.0+\tPN:blastn");
        let fields: Vec<&str> = lines[3].split('\t').collect();
        assert_eq!(fields[0], "Query_1");
        assert_eq!(fields[2], "Subject_1");
        assert_eq!(fields[3], "21");
        // The target is the query: its three letters have no subject letters (I).
        assert_eq!(fields[5], "4H10M3I10M33H");
        let expected: String = (4..27).map(|at| b"AACCGGTTAC"[at % 10] as char).collect();
        assert_eq!(fields[9], expected);
    }

    #[test]
    fn a_minus_hsp_is_reversed_and_complemented() {
        let ops = [GapEditOp::Sub(12)];
        // Query row on the minus strand: the HSP covers the query context offsets 8..20 of 60,
        // so the plus-strand range is 40..52.
        let hsp = hsp((8, 20), (30, 42), -1, 12, &ops);
        let text = document(SamSpec::parse("SQ"), None, hsp, &mut Fixed);
        let lines: Vec<&str> = text.lines().collect();
        let fields: Vec<&str> = lines[3].split('\t').collect();
        assert_eq!(fields[1], "16");
        assert_eq!(fields[3], "41");
        // Subject range [30, 42): reverse clip at the back (30H) and the front (58H).
        assert_eq!(fields[5], "58H12M30H");
        assert_eq!(fields[13], "NM:i:0");
        assert_eq!(fields[14], "PI:f:100.00");
        assert_eq!(fields[9].len(), 12);
    }

    fn query_seq(index: usize, id: &'static str) -> SamSequence<'static> {
        SamSequence {
            record: SamId::Query(QueryRecord {
                id: SeqId::Local(id),
                searched_length: 60,
                offset: 0,
            }),
            key: index as u64,
            index,
            length: 60,
        }
    }

    fn subject_seq(index: usize, id: &'static str) -> SamSequence<'static> {
        SamSequence {
            record: SamId::Subject(SubjectRecord {
                id: SeqId::Local(id),
                searched_length: 100,
                offset: 0,
                use_this_seqid: &[],
            }),
            key: index as u64,
            index,
            length: 100,
        }
    }

    fn run_record() -> RunRecord<'static> {
        RunRecord {
            program: crate::blastinput::app::Program::Blastn,
            db_title: None,
            gapped: true,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1401-1430
    // (`SortHspByMasterStartAscending`): the lowest query coordinate, then the e-value; `list::sort`
    // is stable, so HSPs equal in both keep the order of the Seq-align.
    #[test]
    fn the_default_mode_sorts_a_querys_records_stably_by_start_then_evalue() {
        let ops = [GapEditOp::Sub(10)];
        let mut hsps = Vec::new();
        // (query start, e-value, score): input order 0..4
        for (start, evalue, score) in [
            (30, 1e-5, 11),
            (10, 1e-3, 12),
            (10, 1e-9, 13),
            (30, 1e-5, 14),
            (5, 1.0, 15),
        ] {
            hsps.push(HspRecord {
                evalue,
                score,
                ..hsp((start, start + 10), (0, 10), 1, 10, &ops)
            });
        }
        let mut document = SamDocument::new(SamSpec::default(), "2.17.0+", None);
        document
            .query(
                &run_record(),
                &query_seq(0, "Query_1"),
                &[subject_seq(0, "Subject_1")],
                &hsps,
                &mut NoLetters,
            )
            .unwrap();
        let mut out = Vec::new();
        document.finish(&mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let scores: Vec<&str> = text
            .lines()
            .filter(|line| !line.starts_with('@'))
            .map(|line| line.split('\t').nth(11).unwrap())
            .collect();
        assert_eq!(
            scores,
            ["AS:i:15", "AS:i:13", "AS:i:12", "AS:i:11", "AS:i:14"]
        );
        // `SR` keeps the Seq-align order.
        let mut document = SamDocument::new(SamSpec::parse("SR"), "2.17.0+", None);
        document
            .query(
                &run_record(),
                &query_seq(0, "Query_1"),
                &[subject_seq(0, "Subject_1")],
                &hsps,
                &mut NoLetters,
            )
            .unwrap();
        let mut out = Vec::new();
        document.finish(&mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let scores: Vec<&str> = text
            .lines()
            .filter(|line| !line.starts_with('@'))
            .map(|line| line.split('\t').nth(11).unwrap())
            .collect();
        assert_eq!(
            scores,
            ["AS:i:11", "AS:i:12", "AS:i:13", "AS:i:14", "AS:i:15"]
        );
    }

    // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:445-452
    // (`AddSequence` drops a reference it has): one `@SQ` per reference of the whole run, in order of
    // first appearance, whichever call adds it.
    #[test]
    fn the_sq_lines_list_each_reference_once_across_calls() {
        let ops = [GapEditOp::Sub(10)];
        let hsps = [
            hsp((0, 10), (0, 10), 1, 10, &ops),
            HspRecord {
                subject: 1,
                ..hsp((20, 30), (5, 15), 1, 10, &ops)
            },
        ];
        let subjects = [subject_seq(0, "Subject_1"), subject_seq(1, "Subject_2")];
        for (spec, expected) in [
            (
                SamSpec::default(),
                vec!["@SQ\tSN:Query_1\tLN:60", "@SQ\tSN:Query_2\tLN:60"],
            ),
            (
                SamSpec::parse("SR"),
                vec!["@SQ\tSN:Subject_1\tLN:100", "@SQ\tSN:Subject_2\tLN:100"],
            ),
        ] {
            let mut document = SamDocument::new(spec, "2.17.0+", None);
            for (index, id) in [(0, "Query_1"), (1, "Query_2"), (2, "Query_1")] {
                // the third call names the first query's reference key again
                let mut query = query_seq(index, id);
                query.key = if index == 2 { 0 } else { index as u64 };
                document
                    .query(&run_record(), &query, &subjects, &hsps, &mut NoLetters)
                    .unwrap();
            }
            let mut out = Vec::new();
            document.finish(&mut out).unwrap();
            let text = String::from_utf8(out).unwrap();
            let sq: Vec<&str> = text.lines().filter(|l| l.starts_with("@SQ")).collect();
            assert_eq!(sq, expected);
            assert_eq!(text.lines().filter(|l| !l.starts_with('@')).count(), 6);
            // a second finish writes nothing
            let mut again = Vec::new();
            document.finish(&mut again).unwrap();
            assert!(again.is_empty());
        }
    }

    #[test]
    fn an_empty_document_writes_nothing_and_finish_empties_it() {
        let mut document = SamDocument::new(SamSpec::default(), "2.17.0+", None);
        let mut out = Vec::new();
        document.finish(&mut out).unwrap();
        assert!(out.is_empty());
        assert!(!document.has_records());
    }
}
