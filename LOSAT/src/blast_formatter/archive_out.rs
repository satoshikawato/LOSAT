//! `-outfmt 11` of `blast_formatter` (archive re-emission): NCBI's `CBlastFormat::WriteArchive`
//! and `BlastBuildArchive` (c++/src/algo/blast/format/blast_format.cpp:1066-1133,
//! build_archive.cpp:84-203) over the archive object that `CRemoteBlast::LoadFromArchive` read.
//!
//! Each search object writes one `Blast4-archive` in the encoding named by `ARCHIVE_FORMAT`
//! (`ArchiveEncoding::from_env`: ASN.1 text by default, `xml` for XML, `bin*` for ASN.1 binary,
//! anything else writes nothing). Unlike formats 8, 9 and 12, `-outfmt 11` does not prune the
//! alignments (`PruneSeqalign`) and does not run `BlastFormatter_PreFetchSequenceData`; instead,
//! `-max_target_seqs` and `-num_alignments` update `HitlistSize` in `algorithm-options`, and
//! stored `messages` are never written back into the emitted archive.

use std::cell::RefCell;
use std::collections::BTreeSet;
use std::io::Write;

use super::archive::{invalid_archive, stop_error, Output, SkipWriter, Spill, SPILL_LAG};
use super::refusals::Refusal;
use crate::algorithm::pssm_search::in_pssm::{validation_error, PssmChecks};
use crate::blastinput::app::engine_error;
use crate::ncbi_app::Diag;
use crate::report::archive::fields;
use crate::report::archive::masks::{write_blast4_mask, Frame};
use crate::report::archive::options::{Param, ParamValue, RemoteOptions};
use crate::report::archive::read::options::{BuiltOptions, EProgram};
use crate::report::archive::read::request::{Bioseq, Queries};
use crate::report::archive::read::results::MaskRange;
use crate::report::archive::read::Search;
use crate::report::archive::writer::{
    write_archive_document, write_param, KaBlock, QueryList, RequestSource, ResultsSource,
};
use crate::report::archive::ArchiveEncoding;
use crate::report::seqid::SeqId;
use crate::report::serial::asn_ber::AsnBinaryEncoder;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::write::{write_value, ValueWriter};
use crate::report::serial::asn_types::Value;
use crate::report::serial::encoder::{Encoder, ReplacedChar, SerialError, Traced};
use crate::report::serial::xml::XmlEncoder;
use crate::seqid::accguide::is_local_id;
use crate::seqid::ParsedSeqId;

/// `CPsiBlastValidate::Pssm`'s reading of a `results.pssm` (`PssmChecks`): the PSSM search's
/// archive that `BlastBuildArchive` writes back as `queries pssm`.
fn pssm_checks(pssm_with_params: V<'_>) -> PssmChecks {
    let pssm = pssm_with_params.get("pssm");
    let member = |name: &str| pssm.and_then(|pssm| pssm.get(name));
    let final_data = member("finalData");
    PssmChecks {
        has_pssm: pssm.is_some(),
        missing_scores: final_data
            .and_then(|data| data.get("scores"))
            .is_none_or(|scores| scores.len() == 0),
        missing_freq_ratios: member("intermediateData")
            .and_then(|data| data.get("freqRatios"))
            .is_none_or(|ratios| ratios.len() == 0),
        scaling_factor: final_data
            .and_then(|data| data.get("scalingFactor"))
            .map(|value| value.int()),
        query_is_seq: member("query").map(|query| query.choice().0 == "seq"),
        is_protein: member("isProtein")
            .is_none_or(|value| matches!(value.value, Value::Bool(true))),
    }
}

/// The query part of the emitted `Blast4-queue-search-request`.
enum EmitQueries<'a> {
    Pssm(V<'a>),
    Bioseqs(Vec<Bioseq<'a>>),
    SeqLocs(Vec<&'a ParsedSeqId>),
}

/// The subject part of the emitted `Blast4-queue-search-request`.
enum EmitSubject<'a> {
    Database(&'a [u8]),
    Sequences(Vec<Bioseq<'a>>),
}

/// One `Blast4-mask` of the emitted results (`CRemoteBlast::ConvertToRemoteMasks`).
struct EmitMask<'a> {
    ranges: Vec<MaskRange<'a>>,
    frame: Frame,
}

/// The request of the emitted archive: what `CExportStrategy` makes of the options handle and of
/// the queries and subjects that `CRemoteBlast` read (`writer::write_request_body` writes it).
struct RequestOut<'a> {
    program: &'static str,
    service: &'static str,
    queries: EmitQueries<'a>,
    subject: EmitSubject<'a>,
    algorithm_options: &'a [Param],
    program_options: Vec<Param>,
    step_number: Option<i32>,
}

impl RequestSource for RequestOut<'_> {
    fn program_and_service(&self) -> (&str, &str) {
        (self.program, self.service)
    }

    fn queries(&self) -> QueryList {
        match &self.queries {
            EmitQueries::Pssm(_) => QueryList::Pssm,
            EmitQueries::Bioseqs(bioseqs) => QueryList::Bioseqs(bioseqs.len()),
            EmitQueries::SeqLocs(ids) => QueryList::SeqLocs(ids.len()),
        }
    }

    /// The PSSM and the Bioseqs as read; a Seq-loc `whole` of the Bioseq's first Seq-id.
    fn write_query<E: Encoder>(&self, enc: &mut E, _: QueryList, index: usize) {
        match &self.queries {
            EmitQueries::Pssm(pssm) => write_value(enc, pssm.ty, pssm.value),
            EmitQueries::Bioseqs(bioseqs) => {
                let bioseq = bioseqs[index].value();
                write_value(enc, bioseq.ty, bioseq.value);
            }
            EmitQueries::SeqLocs(ids) => {
                enc.variant("Seq-loc", "whole", 2);
                SeqId::Db(ids[index]).write(enc);
                enc.end_variant();
            }
        }
    }

    fn database(&self) -> Option<&[u8]> {
        match self.subject {
            EmitSubject::Database(name) => Some(name),
            EmitSubject::Sequences(_) => None,
        }
    }

    fn subjects(&self) -> usize {
        match &self.subject {
            EmitSubject::Database(_) => 0,
            EmitSubject::Sequences(bioseqs) => bioseqs.len(),
        }
    }

    fn write_subject<E: Encoder>(&self, enc: &mut E, index: usize) {
        if let EmitSubject::Sequences(bioseqs) = &self.subject {
            let bioseq = bioseqs[index].value();
            write_value(enc, bioseq.ty, bioseq.value);
        }
    }

    fn write_algorithm_options<E: Encoder>(&self, enc: &mut E) {
        for param in self.algorithm_options {
            write_param(enc, param);
        }
    }

    fn has_program_options(&self) -> bool {
        !self.program_options.is_empty()
    }

    fn write_program_options<E: Encoder>(&self, enc: &mut E) {
        for param in &self.program_options {
            write_param(enc, param);
        }
    }

    fn step_number(&self) -> Option<i32> {
        self.step_number
    }
}

/// The results of the emitted archive: the Seq-aligns as read, the masks of
/// `ConvertToRemoteMasks`, the blocks and the statistics that `GetResultSet` read
/// (`writer::write_archive_document` writes them; no `messages`).
struct ResultsOut<'a> {
    alignments: Vec<V<'a>>,
    masks: Vec<EmitMask<'a>>,
    ka_blocks: [KaBlock; 2],
    search_space: i64,
    length_adjustment: i64,
}

impl ResultsSource for &ResultsOut<'_> {
    fn write_alignments<E: Encoder>(
        &mut self,
        enc: &mut E,
        flush: &mut dyn FnMut(&mut E),
    ) -> anyhow::Result<()> {
        // One `ValueWriter` for the Seq-aligns, as `ArchiveDocument::write_aligns`.
        let mut writer = ValueWriter::new();
        for align in &self.alignments {
            enc.element();
            writer.write(enc, align.ty, align.value);
            flush(enc);
        }
        Ok(())
    }

    fn masks(&self) -> usize {
        self.masks.len()
    }

    fn write_mask<E: Encoder>(&self, enc: &mut E, index: usize) {
        let mask = &self.masks[index];
        // The reader's table reads the positions as `TSeqPos` (`Uint4`).
        write_blast4_mask(
            enc,
            mask.frame,
            mask.ranges
                .iter()
                .map(|range| (SeqId::Db(range.id), range.from as u32, range.to as u32)),
        );
    }

    fn ka_blocks(&self) -> [KaBlock; 2] {
        self.ka_blocks
    }

    fn search_stats(&self) -> (i64, i64) {
        (self.search_space, self.length_adjustment)
    }
}

/// Everything needed to serialize one `Blast4-archive` document.
struct ArchiveDoc<'a> {
    request: RequestOut<'a>,
    results: ResultsOut<'a>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1063-1138
/// `CRemoteBlast::ConvertToRemoteMasks` of the per-query mask vectors for `program`: each query's
/// masks in one `Blast4-mask` (`plus1` for a nucleotide query, `notset` for a protein one), a
/// translated query's in one per run of a frame (`FrameNumber2NetworkFrame`).
fn convert_to_remote_masks(
    per_query: Vec<Vec<MaskRange<'_>>>,
    program: EProgram,
) -> Vec<EmitMask<'_>> {
    let translated = program.query_is_translated();
    let default_frame = if program.query_is_nucleotide() || translated {
        Frame::Plus1
    } else {
        Frame::NotSet
    };
    let mut out = Vec::new();
    for query_masks in per_query {
        let Some(first) = query_masks.first() else {
            continue;
        };
        let mut current_frame = first.frame;
        let mut run = Vec::new();
        for range in query_masks {
            if translated && range.frame != current_frame {
                if !run.is_empty() {
                    out.push(EmitMask {
                        ranges: std::mem::take(&mut run),
                        frame: Frame::of_number(current_frame),
                    });
                }
                current_frame = range.frame;
            }
            run.push(range);
        }
        if !run.is_empty() {
            let frame = if translated {
                Frame::of_number(current_frame)
            } else {
                default_frame
            };
            out.push(EmitMask { ranges: run, frame });
        }
    }
    out
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:374-385
// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:477-497
/// `program-options` for a database search: positive `TaxidList` (sorted, deduplicated) or an
/// empty `NegativeTaxidList` (matching `blast_formatter.cpp:381-385`, which constructs the
/// negative list from `taxid_list`).
fn build_program_options(search_options: &BuiltOptions, has_database: bool) -> Vec<Param> {
    if !has_database {
        return Vec::new();
    }
    let mut out = Vec::new();
    if let Some(list) = search_options
        .taxid_list
        .as_ref()
        .filter(|list| !list.is_empty())
    {
        let sorted: BTreeSet<i32> = list.iter().copied().collect();
        out.push(Param::new(
            fields::TAXID_LIST,
            ParamValue::IntegerList(sorted.into_iter().collect()),
        ));
    }
    if search_options
        .negative_taxid_list
        .as_ref()
        .is_some_and(|list| !list.is_empty())
    {
        out.push(Param::new(
            fields::NEGATIVE_TAXID_LIST,
            ParamValue::IntegerList(Vec::new()),
        ));
    }
    out
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-451
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1066-1133
// NCBI reference (598d8ae6): c++/src/algo/blast/format/build_archive.cpp:84-203
/// Formats one archive search object as `-outfmt 11` (`Blast4-archive`).
#[allow(clippy::too_many_arguments)]
pub(super) fn write_archive(
    search: &Search,
    search_options: &BuiltOptions,
    algo_opts: &mut RemoteOptions,
    database: Option<&[u8]>,
    psi_step: Option<Result<i32, ()>>,
    output: &mut Output<'_>,
    buf: &mut Vec<u8>,
    diag: &mut Diag,
) -> anyhow::Result<()> {
    if search_options.program == EProgram::Mapper {
        return Err(Refusal::Program(EProgram::Mapper.task_name()).into());
    }
    let sets = search
        .results
        .alignment_sets(&search.request.queries)
        .map_err(stop_error)?;
    let raw_masks = search
        .results
        .masks(search.request.queries.num_queries(), || {
            search.request.query_is_translated()
        })
        .map_err(stop_error)?;
    let stats = search.results.search_stats(&search.request.service);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2205-2215
    // `GetResultSet` calls `SetMatrixName` for rpsblast and rpstblastn.
    if let Some(matrix) = &stats.matrix {
        algo_opts.set(fields::MATRIX_NAME, ParamValue::Str(matrix.clone()));
    }

    let is_psiblast = search_options.program == EProgram::PsiBlast;
    let step_number = if is_psiblast {
        match psi_step {
            Some(Err(())) => return Err(invalid_archive()),
            Some(Ok(step)) if step != 0 && database.is_some() => Some(step),
            _ => None,
        }
    } else {
        None
    };

    let (program, handle_service) = algo_opts.program_and_service();
    let (queries, subject, service) =
        build_request_endpoints(search, is_psiblast, database, handle_service)?;

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2243-2245
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/build_archive.cpp:95-100,124-130
    // `CBlastAncillaryData` fills the PSI Karlin blocks when `m_Task` is `"psiblast"` and the
    // regular ones otherwise; `s_BuildArchiveAll` reads the PSI blocks when the emitted request's
    // service is `"psi"` (so `psitblastn` reads null pointers and emits `-1, -1, -1`).
    let emitted_is_psi = service == "psi";
    let ka_blocks = if is_psiblast == emitted_is_psi {
        let kb = search.results.karlin_blocks();
        [
            KaBlock {
                lambda: kb.ungapped.0,
                k: kb.ungapped.1,
                h: kb.ungapped.2,
                gapped: false,
            },
            KaBlock {
                lambda: kb.gapped.0,
                k: kb.gapped.1,
                h: kb.gapped.2,
                gapped: true,
            },
        ]
    } else {
        [KaBlock::missing(false), KaBlock::missing(true)]
    };

    let all_alignments = search.results.alignments();
    let alignments: Vec<V<'_>> = sets
        .into_iter()
        .flat_map(|range| all_alignments[range].iter().copied())
        .collect();
    let masks = convert_to_remote_masks(raw_masks, search_options.program);
    let program_options = build_program_options(search_options, database.is_some());

    let doc = ArchiveDoc {
        request: RequestOut {
            program,
            service,
            queries,
            subject,
            algorithm_options: algo_opts.params(),
            program_options,
            step_number,
        },
        results: ResultsOut {
            alignments,
            masks,
            ka_blocks,
            search_space: stats.effective_search_space,
            length_adjustment: stats.length_adjustment,
        },
    };

    let encoding = ArchiveEncoding::from_env();
    if encoding == ArchiveEncoding::Nothing {
        return Ok(());
    }
    // A `FilterString` value that NCBI's `sscanf` did not store, which no later setter set again:
    // the options that NCBI writes hold what its memory held.
    if let Some(what) = algo_opts.undecided() {
        return Err(Refusal::Undecided(what).into());
    }
    output.write_object(diag, |sink, diagnostics| {
        write_archive_streamed(encoding, &doc, buf, sink, diagnostics)
    })
}

/// Builds the `queries`, `subject`, and `service` of the emitted request, validating `results.pssm`
/// and refusing `Queries::SeqLocList` where NCBI resolves locations through its data loaders.
fn build_request_endpoints<'a>(
    search: &'a Search,
    is_psiblast: bool,
    database: Option<&'a [u8]>,
    default_service: &'static str,
) -> anyhow::Result<(EmitQueries<'a>, EmitSubject<'a>, &'static str)> {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:434-446
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:587-598
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:408-420
    if is_psiblast {
        if let Some(pssm) = search.results.reply().get("pssm") {
            if let Some(error) = validation_error(&pssm_checks(pssm)) {
                return Err(engine_error(error));
            }
            let Some(name) = database else {
                return Err(engine_error("Empty reference for CSearchDatabase."));
            };
            if name.is_empty() {
                return Err(engine_error("Error: No database specified"));
            }
            return Ok((EmitQueries::Pssm(pssm), EmitSubject::Database(name), "psi"));
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:609-656
    let queries = match &search.request.queries {
        Queries::SeqLocList(_) => return Err(Refusal::SeqLocQuery.into()),
        Queries::BioseqSet(_) | Queries::Pssm(_) => {
            let bioseqs = search.request.queries.bioseqs();
            let any_local = bioseqs
                .iter()
                .any(|bioseq| bioseq.ids().next().is_some_and(is_local_id));
            if any_local {
                EmitQueries::Bioseqs(bioseqs)
            } else {
                let ids = bioseqs
                    .into_iter()
                    .map(|bioseq| bioseq.ids().next().expect("checked in check_queries"))
                    .collect();
                EmitQueries::SeqLocs(ids)
            }
        }
    };

    let subject = match database {
        Some(name) => {
            if name.is_empty() {
                return Err(engine_error("Error: No database specified"));
            }
            EmitSubject::Database(name)
        }
        None => EmitSubject::Sequences(search.request.subject.bioseqs()),
    };

    Ok((queries, subject, default_service))
}

/// Encodes `doc` into `enc` (`writer::write_archive_document`), spilling finalized bytes after
/// each element of a list when `spill` is present.
fn encode_archive<E: Encoder>(
    enc: &mut E,
    doc: &ArchiveDoc<'_>,
    spill: Option<&RefCell<Spill<'_>>>,
) {
    let mut flush = |enc: &mut E| {
        if let Some(spill) = spill {
            spill.borrow_mut().after_element(enc);
        }
    };
    let written = write_archive_document(enc, &mut flush, &doc.request, &mut &doc.results, &[]);
    written.expect("the read Seq-aligns are written without an error");
}

/// Streams `doc` in `encoding` through the plain encoder with `Spill`; on
/// `SerialError::NonVisibleChar`, re-encodes through `Traced` (no prefix: `blast_formatter` sets
/// none) dropping the bytes already handed on.
fn write_archive_streamed(
    encoding: ArchiveEncoding,
    doc: &ArchiveDoc<'_>,
    buf: &mut Vec<u8>,
    out: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> anyhow::Result<()> {
    buf.clear();
    let (plain, held, written, io_error) = {
        let spill = RefCell::new(Spill::new(&mut *out, SPILL_LAG));
        let plain = run_plain(encoding, doc, buf, Some(&spill));
        let (held, written, io_error) = spill.into_inner().into_parts();
        (plain, held, written, io_error)
    };
    if let Some(error) = io_error {
        return Err(error.into());
    }
    match plain {
        Ok(()) => {
            out.write_all(&held)?;
            out.write_all(buf)?;
            Ok(())
        }
        Err(SerialError::NonVisibleChar(_)) => {
            let mut skip = SkipWriter {
                inner: out,
                skip: written,
            };
            write_archive_posting(encoding, doc, buf, &mut skip, diagnostics)
        }
        Err(error) => Err(unwritable_archive(error)),
    }
}

fn run_plain(
    encoding: ArchiveEncoding,
    doc: &ArchiveDoc<'_>,
    buf: &mut Vec<u8>,
    spill: Option<&RefCell<Spill<'_>>>,
) -> Result<(), SerialError> {
    match encoding {
        ArchiveEncoding::Text => {
            let mut enc = AsnTextEncoder::new(buf);
            encode_archive(&mut enc, doc, spill);
            enc.error().map_or(Ok(()), Err)
        }
        ArchiveEncoding::Binary => {
            let mut enc = AsnBinaryEncoder::new(buf);
            encode_archive(&mut enc, doc, spill);
            enc.error().map_or(Ok(()), Err)
        }
        ArchiveEncoding::Xml => {
            let mut enc = XmlEncoder::archive(buf);
            encode_archive(&mut enc, doc, spill);
            enc.error().map_or(Ok(()), Err)
        }
        ArchiveEncoding::Nothing => Ok(()),
    }
}

#[cold]
#[inline(never)]
fn write_archive_posting(
    encoding: ArchiveEncoding,
    doc: &ArchiveDoc<'_>,
    buf: &mut Vec<u8>,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> anyhow::Result<()> {
    buf.clear();
    let post = |replaced: &ReplacedChar| {
        diagnostics.write_all(&replaced.post(None))?;
        diagnostics.flush()
    };
    match encoding {
        ArchiveEncoding::Text => write_traced_archive(
            Traced::new(AsnTextEncoder::new(buf), 0).posting(post),
            doc,
            writer,
        ),
        ArchiveEncoding::Binary => write_traced_archive(
            Traced::new(AsnBinaryEncoder::new(buf), 0).posting(post),
            doc,
            writer,
        ),
        ArchiveEncoding::Xml => write_traced_archive(
            Traced::xml(XmlEncoder::archive(buf), 0).posting(post),
            doc,
            writer,
        ),
        ArchiveEncoding::Nothing => Ok(()),
    }
}

fn write_traced_archive<E: Encoder>(
    mut enc: E,
    doc: &ArchiveDoc<'_>,
    writer: &mut dyn Write,
) -> anyhow::Result<()> {
    encode_archive(&mut enc, doc, None);
    if let Some(error) = enc.error() {
        return Err(unwritable_archive(error));
    }
    enc.flush_to(writer, true)?;
    Ok(())
}

fn unwritable_archive(error: SerialError) -> anyhow::Error {
    match error {
        SerialError::Unwritten(what) => Refusal::Unwritable(what, 11).into(),
        SerialError::OctetString => Refusal::Unwritable("an OCTET STRING", 11).into(),
        SerialError::NullElement => {
            Refusal::Undecided("an element stated as NULL in a list").into()
        }
        other => other.into(),
    }
}

#[cfg(test)]
#[path = "archive_out_tests.rs"]
mod archive_out_tests;
