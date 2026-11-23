//! TBLASTX's SAM records (outfmt 17), for a TBLASTX search strategy that BLASTN's application runs
//! (the only application whose formatter writes SAM): the report records of each query's final hit
//! list (`algorithm::tblastx::seqalign::query_records_of`) added to the run's [`SamDocument`], one
//! query at a time, and the document written at the end.
//!
//! The Seq-align of a subject is the Std-segs of the search (one per HSP, `BLASTUngappedHspListToSeqAlign`);
//! NCBI's formatter converts it into one Dense-seg, whose record has no score tags, and fails for a
//! subject whose Std-segs change a row's strand (`report::serial::seqalign::translated_segments`):
//! the exception is caught and posted for that Seq-align (below), the subject has no record and no
//! `@SQ` line, and a run whose every Seq-align failed writes no document at all.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:750-755
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/sam.cpp:54-82
//! NCBI reference (598d8ae6): c++/src/objtools/format/cigar_formatter.cpp:84-135

use std::io::Write;

use anyhow::{Context, Result};

use super::blast_engine::TblastxHsp;
use super::report::SubjectLetters;
use super::seqalign::{query_records_of, SeqalignSubjects};
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::sam::{
    database_sam_id, fasta_sam_id, SamDocument, SamId, SamLetters, SamRow, SamSequence, SamSpec,
};

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:815-830
// NCBI BLAST+ 2.17.0 reports it at line 827 (the oracle's probes of F2 fix round 3, `F/f2/tmp/n1/
// tblastx.sub.ss_.err`).
/// Where `CreateDensegFromStdseg` throws when the Std-segs of a Seq-align put a row on different
/// strands.
const INCONSISTENT_STRANDS: ThrowSite = ThrowSite {
    file: "src/objects/seq/../seqalign/Seq_align.cpp",
    line: 827,
    kind: "CSeqalignException::eInvalidInputAlignment",
    function: None,
};

// NCBI reference (598d8ae6): c++/src/objtools/format/cigar_formatter.cpp:110-125
// NCBI reference (598d8ae6): c++/src/objtools/format/cigar_formatter.cpp:28 (the error code of the
// file, `Objtools_Fmt_CIGAR`: 1439)
/// The post of the formatter's catch-all handler (error code 1 of the file) for a Seq-align whose
/// Std-segs change a strand, in the application `app` (`Error: (1439.1) [<app>] [...] Exception:
/// <the exception's report>`), in LOSAT's words.
pub(crate) fn inconsistent_strands_post(app: &str) -> Vec<u8> {
    let exception = NcbiException::new(
        INCONSISTENT_STRANDS,
        "CreateDensegFromStdseg(): Inconsistent strands!",
    );
    crate::ncbi_exception::error_post_line(
        Some("1439.1"),
        Some(app),
        &[
            b"[CCIGAR_Formatter::x_FormatAlignmentRows] Exception: ".as_slice(),
            &exception.ncbi_what(),
        ]
        .concat(),
    )
}

/// The letters of the queries and the subjects for SEQ (`SQ`): the searched letters, in the
/// coordinates of the records.
struct TblastxSamLetters<'a> {
    /// The searched query records (with `-query_loc`, the interval's letters).
    queries: &'a [FastaRecord],
    query_placements: &'a Placements,
    /// The subjects as the Seq-aligns name them: where each lies in its record (`-subject_loc`).
    subjects: &'a SeqalignSubjects<'a>,
    subject_letters: SubjectLetters<'a>,
}

impl SamLetters for TblastxSamLetters<'_> {
    fn letters(
        &mut self,
        row: SamRow,
        index: usize,
        from: usize,
        to: usize,
        out: &mut Vec<u8>,
    ) -> Result<()> {
        match row {
            SamRow::Query => {
                let offset = self.query_placements.offset(index);
                let record = self.queries[index].seq();
                let letters = from
                    .checked_sub(offset)
                    .zip(to.checked_sub(offset))
                    .and_then(|(from, to)| record.get(from..to))
                    .context("a SAM record lies outside the searched query letters")?;
                out.extend(letters.iter().map(u8::to_ascii_uppercase));
            }
            SamRow::Subject => {
                let offset = match self.subjects {
                    SeqalignSubjects::Records { placements, .. } => placements.offset(index),
                    SeqalignSubjects::Database { .. } => 0,
                };
                let length = self.subject_letters.length(index);
                let (from, to) = from
                    .checked_sub(offset)
                    .zip(to.checked_sub(offset))
                    .filter(|&(from, to)| from <= to && to <= length)
                    .context("a SAM record lies outside the searched subject letters")?;
                self.subject_letters.letters_of(index, from, to, out);
            }
        }
        Ok(())
    }
}

/// Writes the SAM document of the queries (`queries`, their final hits `hits` in report order) to
/// `writer` at the end, each query's warnings before its records and NCBI's post of each Seq-align
/// that its formatter fails on after them (on the warnings' sink). `spec` is the custom
/// specification of `-outfmt 17`, `command_line` the `CL` tag of the `@PG` line (`None`: none),
/// `app` the running application's name (the prefix of the posts), `queries_before` the queries that
/// earlier calls wrote (the numbering of the references).
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_tblastx_sam(
    spec: SamSpec,
    command_line: Option<Vec<u8>>,
    writer: &mut dyn Write,
    hits: &[TblastxHsp],
    queries: &[FastaRecord],
    query_placements: &Placements,
    queries_before: usize,
    subjects: &SeqalignSubjects<'_>,
    subject_letters: SubjectLetters<'_>,
    warnings: &mut QueryWarnings<'_>,
    app: &str,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Tblastx,
        rps: false,
        db_title: match subjects {
            SeqalignSubjects::Records { .. } => None,
            SeqalignSubjects::Database { title, .. } => Some(title),
        },
        gapped: false,
    };
    let mut document = SamDocument::new(spec, crate::wording::LOSAT_VERSION, command_line)
        .for_program(Program::Tblastx);
    let mut letters = TblastxSamLetters {
        queries,
        query_placements,
        subjects,
        subject_letters,
    };
    let mut subject_records: Vec<SubjectRecord<'_>> = Vec::new();
    let mut hsp_records: Vec<HspRecord<'_>> = Vec::new();
    // The records do not read the engine's counts (a TBLASTX Seq-align has no score tags).
    let mut counts: Vec<(u32, u32)> = Vec::new();
    let mut next = 0usize;
    for (q_idx, query) in queries.iter().enumerate() {
        // The query's warnings come before its records (`QueryWarnings`).
        warnings.before_query(q_idx, &mut *writer)?;
        let start = next;
        while next < hits.len() && hits[next].hit.q_idx as usize == q_idx {
            next += 1;
        }
        if start == next {
            // A query without alignments writes nothing (no unmapped record).
            continue;
        }
        if counts.len() < next - start {
            counts.resize(next - start, (0, 0));
        }
        query_records_of(
            &hits[start..next],
            &counts[..next - start],
            subjects,
            subject_letters,
            &mut subject_records,
            &mut hsp_records,
        )?;
        // The subject of each record: its index in the call's subjects, in the order of the runs.
        let mut subject_indexes: Vec<usize> = Vec::with_capacity(subject_records.len());
        let mut last = None;
        for hsp in &hits[start..next] {
            if last != Some(hsp.hit.s_idx) {
                last = Some(hsp.hit.s_idx);
                subject_indexes.push(hsp.hit.s_idx as usize);
            }
        }
        let mut sam_subjects: Vec<SamSequence<'_>> = Vec::with_capacity(subject_records.len());
        for (record, &s_idx) in subject_records.iter().zip(&subject_indexes) {
            let (id, key, length) = match subjects {
                SeqalignSubjects::Records {
                    records,
                    placements,
                } => (
                    fasta_sam_id(&records[s_idx]),
                    s_idx as u64,
                    placements.length(s_idx, record.searched_length as usize) as i32,
                ),
                SeqalignSubjects::Database {
                    ordinals,
                    oids,
                    parsed,
                    ..
                } => (
                    database_sam_id(
                        parsed.get(s_idx).and_then(Option::as_deref),
                        *oids,
                        ordinals[s_idx],
                    ),
                    u64::from(ordinals[s_idx]),
                    record.searched_length,
                ),
            };
            sam_subjects.push(SamSequence {
                record: SamId::Subject(SubjectRecord { id, ..*record }),
                key,
                index: s_idx,
                length,
            });
        }
        let query_record = QueryRecord {
            id: fasta_sam_id(query),
            searched_length: i32::try_from(query.seq().len())?,
            offset: i32::try_from(query_placements.offset(q_idx))?,
        };
        let query_sequence = SamSequence {
            record: SamId::Query(query_record),
            key: (queries_before + q_idx) as u64,
            index: q_idx,
            length: query_placements.length(q_idx, query.seq().len()) as i32,
        };
        let failed = document.query(
            &run,
            &query_sequence,
            &sam_subjects,
            &hsp_records,
            &mut letters,
        )?;
        // NCBI's formatter posts each Seq-align that it cannot convert while it formats the query.
        if failed > 0 {
            writer.flush()?;
            let post = inconsistent_strands_post(app);
            for _ in 0..failed {
                warnings.sink.write_all(&post)?;
            }
            warnings.sink.flush()?;
        }
    }
    debug_assert_eq!(next, hits.len(), "the final list is in query order");
    document.finish(writer)?;
    Ok(())
}

#[cfg(test)]
#[path = "sam_tests.rs"]
mod tests;
