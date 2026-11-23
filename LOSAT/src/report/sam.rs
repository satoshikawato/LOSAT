//! SAM (outfmt 17) of the Seq-aligns of a nucleotide search: NCBI's `CBlast_SAM_Formatter`
//! (`CSAM_Formatter`, `CSAM_CIGAR_Formatter`, `CCIGAR_Formatter`), which keeps every record until
//! the formatter is destroyed and writes the header last. One [`SamDocument`] holds a run's
//! records as rendered bytes, its `@SQ` lines and the references it has listed; [`SamDocument::finish`]
//! writes the document.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:313-316,405-441, :405-442
//!
//! The record of a Seq-align is the one NCBI's `CSAM_CIGAR_Formatter::AddRow` makes from its
//! Dense-seg (`CCIGAR_Formatter::x_FormatLine` for the CIGAR). A gapped search has a Seq-align per
//! HSP, whose Dense-seg is the one that outfmt 8 and 12 write
//! (`report::serial::seqalign::gapped_segments`); an ungapped search has a Seq-align per subject,
//! whose Dense-diags (one per HSP) the formatter joins into one Dense-seg
//! (`report::serial::seqalign::ungapped_segments`), so a subject has one record. A translated
//! ungapped search (TBLASTX's strategy run by BLASTN's application) has a Seq-align of Std-segs per
//! subject, which the formatter converts into one Dense-seg (`translated_segments`): a subject whose
//! Std-segs change a row's strand has no record (the conversion throws; `SamDocument::query` counts
//! it), and a record has no score tags (the Seq-align holds no scores).

use std::collections::HashSet;
use std::io::{self, Write};

use anyhow::Result;

use super::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use super::seqid::SeqId;
use super::serial::seqalign::{
    gapped_segments, seq_align_hsps, seqalign_evalue, translated_segments, ungapped_segments,
    untranslated, Segments, Strand, GAP_VALUE,
};

/// The words of the custom specification of outfmt 17 that NCBI reads.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/sam.cpp:84-104
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

/// The Seq-id that SAM names a FASTA record (a query, or a `-subject` record) by: the local id, or
/// with `-parse_deflines` the id of the record's Bioseq that `AdjustSeqIdType` with
/// `eGetId_ForceAcc` finds (`db_subject::forced_accession`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:122-132
pub fn fasta_sam_id(record: &crate::blastinput::fasta_reader::FastaRecord) -> SeqId<'_> {
    match record.believed() {
        Some(ids) => SeqId::Db(crate::report::db_subject::forced_accession(
            &ids.best, &ids.ids,
        )),
        None => SeqId::Local(&record.local_id),
    }
}

/// The Seq-id that SAM names a database sequence by: an id of its Bioseq where it has parsed
/// Seq-ids (`ParsedSubject::sam_id`), else its `BL_ORD_ID` (the OID of `ordinal`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:122-132
pub fn database_sam_id<'a>(
    subject: Option<&'a crate::report::db_subject::ParsedSubject>,
    oids: crate::seqdb::OidMap<'_>,
    ordinal: u32,
) -> SeqId<'a> {
    match subject {
        Some(subject) => SeqId::Db(subject.sam_id()),
        None => SeqId::BlOrdId(oids.oid(ordinal)),
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
#[cfg(test)]
pub struct NoLetters;

#[cfg(test)]
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
    /// The lesser of `GetSeqStart` and `GetSeqStop` of the query row (row 0), then the Seq-align's e-value: its
    /// last HSP's (`SortHspByMasterStartAscending`).
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
    /// The CIGAR's operations and counts (`TSeqPos`).
    ops: Vec<(u8, u32)>,
}

/// The document of a run.
pub struct SamDocument {
    spec: SamSpec,
    /// The `@PG` line without its CL tag and newline: `@PG\tID:0\tVN:<version>`.
    program: Vec<u8>,
    /// The `PN` tag of the `@PG` line: the options' program (`blastn`, or `tblastx` for a TBLASTX
    /// strategy that BLASTN's application runs).
    program_name: &'static str,
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
    pub fn new(spec: SamSpec, version: &str, command_line: Option<Vec<u8>>) -> Self {
        Self {
            spec,
            program: format!("@PG\tID:0\tVN:{version}").into_bytes(),
            program_name: "blastn",
            command_line: command_line.filter(|line| !line.is_empty()),
            sq: Vec::new(),
            seen: HashSet::new(),
            body: ChunkedBytes::default(),
            scratch: Scratch::default(),
        }
    }

    /// The document of the options' program `program` (its `PN` tag): BLASTN's unless the strategy
    /// that BLASTN's application runs is another program's.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2460-2466
    pub fn for_program(mut self, program: crate::blastinput::app::Program) -> Self {
        self.program_name = program.ncbi_name();
        self
    }

    /// Adds the records of one query's final hit list: `hsps`, grouped by subject in the
    /// Seq-align order, `subjects` indexed by `HspRecord::subject`. `run` says whether the
    /// subjects are remapped (`-subject`) and whether the search is gapped (a Seq-align and a
    /// record per HSP) or ungapped (a Seq-align and a record per subject: its run of HSPs).
    /// Returns the number of Seq-aligns that the formatter could not convert (a translated search's
    /// subject whose Std-segs change a row's strand: no record, one post each, `0` for the other
    /// searches).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/sam.cpp:54-82
    pub fn query(
        &mut self,
        run: &RunRecord<'_>,
        query: &SamSequence<'_>,
        subjects: &[SamSequence<'_>],
        hsps: &[HspRecord<'_>],
        letters: &mut dyn SamLetters,
    ) -> Result<usize> {
        if hsps.is_empty() {
            return Ok(0);
        }
        let mut scratch = std::mem::take(&mut self.scratch);
        scratch.rows.clear();
        scratch.buf.clear();
        let result = (|| -> Result<usize> {
            let mut failed = 0;
            for align in seq_align_hsps(run, hsps) {
                let subject = align[0].subject as usize;
                if !self.render(run, query, &subjects[subject], align, letters, &mut scratch)? {
                    failed += 1;
                    continue;
                }
                scratch.rows.last_mut().expect("a row was rendered").subject = subject;
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
            Ok(failed)
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

    /// Renders the record of one Seq-align into `scratch` (`CSAM_CIGAR_Formatter::AddRow` after
    /// `CCIGAR_Formatter::x_FormatLine`): `align` is a gapped search's HSP, or an ungapped
    /// search's HSPs of one subject in the Seq-align's order. `false`: the formatter cannot convert
    /// the Seq-align (a translated search's Std-segs on different strands of a row), nothing is
    /// rendered.
    fn render(
        &self,
        run: &RunRecord<'_>,
        query: &SamSequence<'_>,
        subject: &SamSequence<'_>,
        align: &[HspRecord<'_>],
        letters: &mut dyn SamLetters,
        scratch: &mut Scratch,
    ) -> Result<bool> {
        let (SamId::Query(query_record), SamId::Subject(subject_record)) =
            (&query.record, &subject.record)
        else {
            unreachable!("the first row is the query and the second the subject")
        };
        // A translated search's Seq-align is made of Std-segs, and its Dense-seg of `CreateDensegFromStdseg`.
        let translated = !run.gapped && !untranslated(run.program);
        if run.gapped {
            debug_assert_eq!(align.len(), 1, "a gapped search's Seq-align holds one HSP");
            gapped_segments(
                run,
                query_record,
                subject_record,
                &align[0],
                &mut scratch.segments,
            );
        } else if translated {
            if translated_segments(query_record, subject_record, align, &mut scratch.segments)
                .is_err()
            {
                return Ok(false);
            }
        } else {
            ungapped_segments(query_record, subject_record, align, &mut scratch.segments);
        }
        // The tags read the first of each score of the Seq-align (`AddRow`): its first HSP's. For
        // an ungapped search `CreateDensegFromDendiag` joins the HSPs' score lists into the
        // Dense-seg's scores and `Print` swaps them into the Seq-align's, each list with a score,
        // an e-value, a bit score and a number of identities (a BLASTN HSP has no sum statistics).
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/sam.cpp:71-74
        // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:247-279
        let hsp = &align[0];
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
        // The reference row and the target row (`FormatByTargetRow` with `m_refRow`).
        let subject_ref = self.spec.subject_as_reference;
        let (tgt_row, ref_row) = if subject_ref { (0, 1) } else { (1, 0) };
        // `StrandSign` of `row`: the row's strand in the first segment.
        let reverse = (strands[tgt_row] == Strand::Minus) != (strands[ref_row] == Strand::Minus);
        let ref_minus = strands[ref_row] == Strand::Minus;
        let numseg = lens.len();
        // The CIGAR (`CCIGAR_Formatter::x_FormatLine`, widths 1), merging neighbours of the
        // same type; the ranges of the rows are the unions of their segments. The counts and the
        // aligned length (`GetAlignLength` with false, the lengths of the segments without a gap) are
        // `TSeqPos` sums, which wrap at 2^32: an ungapped subject's record adds up all its HSPs,
        // whose diagonals may cover the same letters many times (a long tandem repeat).
        // NCBI reference (598d8ae6): c++/src/objtools/format/cigar_formatter.cpp:164, :324-346
        // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:1867-1892
        let ops = &mut scratch.ops;
        ops.clear();
        let mut ref_range = (i64::MAX, i64::MIN);
        let mut tgt_range = (i64::MAX, i64::MIN);
        let mut aligned: u32 = 0;
        for idx in 0..numseg {
            let seg = if ref_minus { numseg - idx - 1 } else { idx };
            let len = i64::from(lens[seg]);
            let count = lens[seg] as u32;
            let ref_start = start_of(seg, ref_row);
            let tgt_start = start_of(seg, tgt_row);
            let ref_seq = ref_start != GAP_VALUE;
            let tgt_seq = tgt_start != GAP_VALUE;
            let kind = if tgt_seq && !ref_seq {
                b'I'
            } else if !tgt_seq && ref_seq {
                b'D'
            } else if tgt_seq && ref_seq {
                aligned = aligned.wrapping_add(count);
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
                Some((last, sum)) if *last == kind => *sum = sum.wrapping_add(count),
                _ => ops.push((kind, count)),
            }
        }
        // NM: `m_NumDif`, an `int`, adds the count of each operation other than M (`AddSegment`).
        // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:150-157
        let num_dif = ops
            .iter()
            .filter(|(kind, _)| *kind != b'M')
            .fold(0i32, |sum, &(_, count)| sum.wrapping_add(count as i32));
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
            push_int(buf, i64::from(count));
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
        // The tags: AS, EV, NM, PI, BS, in the order of the code (`AS + EV + NM + PI + BS`). The
        // Seq-align that `AddRow` reads holds no score where it is a translated search's (its
        // scores are the Std-segs'): NM alone.
        // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:247-253
        if !translated {
            buf.extend_from_slice(b"\tAS:i:");
            push_int(buf, i64::from(hsp.score));
            buf.extend_from_slice(b"\tEV:f:");
            push_g(buf, seqalign_evalue(hsp.evalue));
        }
        buf.extend_from_slice(b"\tNM:i:");
        push_int(buf, i64::from(num_dif));
        if !translated {
            buf.extend_from_slice(b"\tPI:f:");
            // `len`, an `int`, is `GetAlignLength` of `aln` with false: from 2^31 on the length is negative, and so is PI
            // (NCBI BLAST+ 2.17.0 writes `PI:f:-0.01` for a subject whose HSPs add up to 2152959968).
            // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:267-276
            let len = aligned as i32;
            let identity = if hsp.num_ident == len {
                100.0
            } else {
                f64::min(99.99, 100.0 * f64::from(hsp.num_ident) / f64::from(len))
            };
            push_fixed2(buf, identity);
            buf.extend_from_slice(b"\tBS:f:");
            push_g(buf, hsp.bit_score);
        }
        buf.push(b'\n');
        let end = buf.len();
        // The sort reads the last e-value of the Seq-align's scores (`GetAlnScores`), which a
        // translated search's Seq-align (Std-segs, sorted before the conversion) holds at its
        // first Std-seg; the start of the sorted row is the lowest start of the Std-segs' query rows
        // (`GetSeqStart` and `GetSeqStop` of a Std-seg alignment combine the rows' ranges).
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:229-270
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:815-845
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1402-1430
        // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:152-300
        let (start, scored) = if translated {
            (translated_row_start(segments, offsets[0]), &align[0])
        } else {
            (
                query_row_start(segments, offsets[0]),
                align.last().expect("a Seq-align holds an HSP"),
            )
        };
        scratch.rows.push(Rendered {
            start,
            evalue: seqalign_evalue(scored.evalue),
            bytes: begin..end,
            subject: 0,
        });
        Ok(true)
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
        out.write_all(b"\tPN:")?;
        out.write_all(self.program_name.as_bytes())?;
        out.write_all(b"\n")?;
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

/// The lesser of `GetSeqStart` and `GetSeqStop` of the Seq-align's query row (row 0), in record coordinates
/// (`offset` added to the starts). Both read the row's strand in the first segment: on the plus
/// strand the start is the first segment's that is not a gap and the stop the end of the last
/// one, on the minus strand the other way round. A gapped HSP's row runs one way, so this is its
/// lowest coordinate; an ungapped search's segments are its subject's HSPs in e-value order, and
/// a lower HSP between the first and the last one is not read. The positions are `TSeqPos`
/// (unsigned).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:165-238
fn query_row_start(segments: &Segments, offset: i32) -> i32 {
    let starts = segments.starts();
    let lens = segments.lens();
    let mut rows = (0..lens.len()).filter(|&seg| starts[2 * seg] >= 0);
    let Some(first) = rows.next() else {
        return 0;
    };
    let last = rows.next_back().unwrap_or(first);
    let (start_seg, stop_seg) = if segments.strands()[0] == Strand::Minus {
        (last, first)
    } else {
        (first, last)
    };
    let start = (starts[2 * start_seg] + offset) as u32;
    let stop = ((starts[2 * stop_seg] + offset + lens[stop_seg]) as u32).wrapping_sub(1);
    start.min(stop) as i32
}

/// The lowest start of the query row (row 0) of a translated search's Seq-align, in record
/// coordinates: its Std-segs' intervals (`segments` holds one segment each) start there, and
/// `GetSeqStart` and `GetSeqStop` of the Seq-align are the lowest start and the highest end of the
/// intervals; the lesser of them is the lowest start. A `TSeqPos` (unsigned).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:152-300
fn translated_row_start(segments: &Segments, offset: i32) -> i32 {
    let starts = segments.starts();
    (0..segments.lens().len())
        .map(|seg| (starts[2 * seg] + offset) as u32)
        .min()
        .unwrap_or(0) as i32
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

/// `%.2f` (`NStr::DoubleToString` of `value` with 2 digits).
fn push_fixed2(out: &mut Vec<u8>, value: f64) {
    out.extend_from_slice(format!("{value:.2}").as_bytes());
}

/// C's `%g` (precision 6), which `NStr::DoubleToString` of `value` prints.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2316-2340
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
#[path = "sam_tests.rs"]
mod tests;
