//! The `Blast4-archive` document, written element by element to the sink.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/build_archive.cpp:84-152
//! (`s_BuildArchiveAll`: the request of the search strategy with the version as `ident`,
//! the results with the alignments of every query of the batch, the masks, the Karlin
//! blocks and statistics of the batch's first query, then the messages).
//!
//! ```text
//! Blast4-archive ::= {
//!   request { ident "2.17.0+", body queue-search { program, service, queries bioseq-set,
//!             subject, algorithm-options, program-options } },
//!   results { alignments, masks, ka-blocks, search-stats },
//!   messages { ... } }
//! ```
//! LOSAT holds no more than NCBI does for a batch: the encoder's buffer holds the line that is
//! being written and everything else has gone to the sink (`Encoder::flush_to`).

use std::io::Write;

use anyhow::Result;

use super::bioseq::{write_bioseq, write_seq_loc, BioseqInput};
use super::fields;
use super::masks::{write_mask, MaskEntry};
use super::options::{Param, ParamValue};
use super::{ArchiveEncoding, ArchiveHandle};
use crate::report::diag_log::Message;
use crate::report::serial::asn_ber::AsnBinaryEncoder;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::encoder::{Encoder, ReplacedChar, Traced};
use crate::report::serial::pssm::PssmWithParameters;
use crate::report::serial::xml::XmlEncoder;

/// The bytes that an encoder holds before its buffer is handed to the sink.
const FLUSH_BYTES: usize = 32 * 1024;

/// A Karlin-Altschul block of the results (`Blast4-ka-block`); `-1` three times is the block
/// that does not exist.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KaBlock {
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
    pub gapped: bool,
}

impl KaBlock {
    /// The block `s_Convert_to_CBlast_ka_block` writes for a missing `Blast_KarlinBlk`.
    pub fn missing(gapped: bool) -> Self {
        Self {
            lambda: -1.0,
            k: -1.0,
            h: -1.0,
            gapped,
        }
    }
}

/// The queries of a request (`Blast4-queries`): the batch's Bioseqs (a FASTA query's), or the
/// PSSM of a PSSM query (`-in_pssm`) as NCBI's application read it.
#[derive(Clone, Copy)]
pub enum RequestQueries<'a> {
    /// The batch's queries, in order.
    Bioseqs(&'a [BioseqInput<'a>]),
    /// The `PssmWithParameters` of psiblast's or tblastn's `-in_pssm`.
    Pssm(&'a PssmWithParameters),
}

/// What the request of one query batch reads besides the run's options (`ArchiveRequest`): the
/// archive's and the search strategy's (`requests::RunRequests` makes it from the batch's
/// records).
#[derive(Clone, Copy)]
pub struct RequestInputs<'a> {
    /// The batch's queries.
    pub queries: RequestQueries<'a>,
    /// `-subject`: every subject (the request's `database` is `None` then).
    pub subjects: &'a [BioseqInput<'a>],
    /// `RequiredStart` and `RequiredEnd` (0-based) of a `-query_loc` that is not the first
    /// query's whole length.
    pub required: Option<(i32, i32)>,
    /// The `LCaseMask` parameters, in order.
    pub lcase_masks: &'a [MaskEntry],
}

impl RequestInputs<'_> {
    /// Whether `encoding` writes the strings of the request's Bioseqs and Seq-locs as they are
    /// (`BioseqInput::written_as_is`); the `LCaseMask` parameters, the masks and the Seq-aligns
    /// of the archive hold the same Seq-ids (the subjects of `-db` have the database's, which
    /// its reader refuses or replaces where NCBI would post). The strings of a PSSM query are
    /// those of the PSSM's own text document (`pssm::asn_text_with_karlin`: one owner of what the
    /// ASN.1 text stream replaces): an `-in_pssm` file's are visible ASCII (NCBI's reader replaces
    /// any other byte before the request is made), a title that `-in_msa` takes from a defline
    /// need not be. A PSSM with a list element stated as `NULL` (`-in_pssm`) is written through
    /// `Traced` in every stream, which throws there as NCBI's does (`encoder::NullElement`).
    pub fn written_as_is(&self, encoding: ArchiveEncoding) -> bool {
        let (queries, pssm_as_is) = match self.queries {
            RequestQueries::Bioseqs(queries) => (queries, true),
            RequestQueries::Pssm(pssm) => (
                &[][..],
                !super::pssm_holds_null_element(pssm)
                    && match encoding {
                        ArchiveEncoding::Text => {
                            crate::report::serial::pssm::asn_text_with_karlin(pssm, None).is_ok()
                        }
                        // A UTF8String's control byte (NCBI's text reader replaced those of the
                        // VisibleStrings), which NCBI's XML stream replaces and posts; in both
                        // streams a delayed member of another format whose read posts a replaced
                        // byte (`Encoder::delayed_posts`, written through `Traced`).
                        ArchiveEncoding::Xml | ArchiveEncoding::Binary => {
                            super::pssm_written_as_is(encoding, pssm)
                        }
                        ArchiveEncoding::Nothing => true,
                    },
            ),
        };
        pssm_as_is
            && queries
                .iter()
                .chain(self.subjects)
                .all(|bioseq| bioseq.written_as_is(encoding))
    }
}

/// What one archive holds, besides the request's options (`ArchiveHandle::request`).
pub struct BatchArchive<'a> {
    /// The batch's request.
    pub request: RequestInputs<'a>,
    /// The `masks` of the results.
    pub masks: &'a [MaskEntry],
    /// The ungapped and the gapped block of the batch's first query.
    pub ka_blocks: [KaBlock; 2],
    /// `search-stats`: the effective search space (written twice) and the length adjustment.
    pub search_space: i64,
    pub length_adjustment: i64,
    /// The batch's context 0 when the batch was searched, which its `EffectiveSearchSpace`
    /// option reads (`RequestLength::space`).
    pub searched_space: Option<SearchedContext0>,
    /// psiblast's round (`ArchiveRound`); `None` for the other programs.
    pub round: Option<ArchiveRound>,
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:259-271
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1066-1133
/// What psiblast's archive of one round sets in its request beyond the run's (`WriteArchive` with
/// the round's number, after the round's search): the `Web_StepNumber` format option is the round's
/// number with `-db` (`BlastBuildArchive`'s `num_iters`; the bl2seq archive of `-subject` has
/// none), and a PSSM search's set-up leaves its own composition-based mode in the options
/// (`PsiBlastSetupScoreBlock` resets 2 and 3 to 1 for the query's later rounds;
/// `algorithm::psiblast::in_pssm::position_based_mode`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ArchiveRound {
    /// `itr.GetIterationNumber()`, from 1.
    pub number: i32,
    /// The `CompositionBasedStats` option of a PSSM search (`None`: the command line's).
    pub comp_based_stats: Option<i32>,
}

/// What the `EffectiveSearchSpace` option of a searched batch reads of the batch's context 0
/// (`RequestLength::space`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SearchedContext0 {
    /// The search's effective search space of the context (`-searchsp` in an invalid one too).
    pub space: i64,
    /// The context's query length, 0 for an invalid context.
    pub query_length: i64,
}

/// The encoder, the sink and the first I/O error: the documents are written by calls that
/// do not return one.
pub struct ArchiveOut<'s, E: Encoder> {
    pub enc: E,
    sink: &'s mut dyn Write,
    io_error: Option<std::io::Error>,
}

impl<'s, E: Encoder> ArchiveOut<'s, E> {
    /// The encoder `enc` over `sink`.
    fn new(enc: E, sink: &'s mut dyn Write) -> Self {
        ArchiveOut {
            enc,
            sink,
            io_error: None,
        }
    }

    /// Hands the encoder's final bytes to the sink once enough have gathered. Call it after
    /// each element of a long list.
    pub fn flush_if_large(&mut self) {
        if self.io_error.is_none() && self.enc.buffered() >= FLUSH_BYTES {
            if let Err(error) = self.enc.flush_to(self.sink, false) {
                self.io_error = Some(error);
            }
        }
    }

    /// Runs `write` with the encoder and a function that hands the encoder's final bytes to
    /// the sink (what `flush_if_large` does), for a writer that holds the encoder while it
    /// writes an element.
    pub fn write_with<R>(&mut self, write: impl FnOnce(&mut E, &mut dyn FnMut(&mut E)) -> R) -> R {
        let ArchiveOut {
            enc,
            sink,
            io_error,
        } = self;
        let mut flush = |enc: &mut E| {
            if io_error.is_none() && enc.buffered() >= FLUSH_BYTES {
                if let Err(error) = enc.flush_to(&mut **sink, false) {
                    *io_error = Some(error);
                }
            }
        };
        write(enc, &mut flush)
    }

    pub(super) fn finish(mut self) -> Result<()> {
        if self.io_error.is_none() {
            if let Err(error) = self.enc.flush_to(self.sink, true) {
                self.io_error = Some(error);
            }
        }
        if let Some(error) = self.enc.error() {
            return Err(error.into());
        }
        match self.io_error {
            Some(error) => Err(error.into()),
            None => Ok(()),
        }
    }
}

impl<E: Encoder, P: FnMut(&ReplacedChar) -> std::io::Result<()>> ArchiveOut<'_, Traced<E, P>> {
    /// `finish` of a document whose stream threw (`encoder::NullElement`, the application
    /// `application`): the sink gets the bytes that NCBI's stream flushed (with a post function,
    /// each replaced byte posted after the bytes NCBI had written before it), `collect` the bytes
    /// replaced before the throw that no post function posted, then NCBI's error ends the request
    /// (`cli::NativeError::null_element`).
    fn finish_thrown(
        mut self,
        collect: impl FnOnce(&[ReplacedChar]) -> std::io::Result<()>,
        application: &str,
        thrown: &crate::report::serial::encoder::NullElement,
    ) -> Result<()> {
        if self.io_error.is_none() {
            if let Err(error) = self.enc.flush_to(self.sink, true) {
                self.io_error = Some(error);
            }
        }
        let replaced = self.enc.take_replaced_chars();
        if !replaced.is_empty() {
            collect(&replaced)?;
        }
        match self.io_error {
            Some(error) => Err(error.into()),
            None => Err(crate::cli::NativeError::null_element(application, thrown).into()),
        }
    }

    /// `finish` of an archive written through `Traced` with a post function, or `finish_thrown`
    /// where its stream threw: then the error is NCBI's (`cli::NativeError`), which no other
    /// error of `finish` is.
    fn finish_archive(mut self, application: &str) -> Result<()> {
        match self.enc.take_thrown() {
            Some(null_element) => self.finish_thrown(|_| Ok(()), application, &null_element),
            None => self.finish(),
        }
    }
}

/// The Seq-aligns of a batch: written as the elements of the `alignments` list.
pub trait AlignSource {
    fn write<E: Encoder>(&mut self, out: &mut ArchiveOut<'_, E>) -> Result<()>;
}

fn write_param<E: Encoder>(enc: &mut E, param: &Param) {
    enc.element();
    enc.begin_sequence("Blast4-parameter");
    enc.member("name", 0);
    enc.visible_str(param.name().as_bytes());
    enc.end_member();
    enc.member("value", 1);
    let value_type = param.value.value_type();
    enc.variant("Blast4-value", value_type.variant(), value_type.tag());
    match &param.value {
        ParamValue::Evalue(value) => {
            enc.variant("Blast4-cutoff", "e-value", 0);
            enc.real(*value);
            enc.end_variant();
        }
        ParamValue::Integer(value) => enc.int(*value),
        ParamValue::BigInteger(value) => enc.big_int(*value),
        ParamValue::Boolean(value) => enc.boolean(*value),
        ParamValue::Real(value) => enc.real(*value),
        ParamValue::Str(value) => enc.visible_str(value.as_bytes()),
        ParamValue::IntegerList(values) => {
            enc.begin_list(false);
            for &value in values {
                enc.element();
                enc.int(value);
            }
            enc.end_list();
        }
        // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:493-497
        ParamValue::StrandType(value) => enc.enumerated(
            match value {
                1 => "forward-strand",
                2 => "reverse-strand",
                _ => "both-strands",
            },
            *value,
        ),
    }
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
}

fn write_messages<E: Encoder>(enc: &mut E, messages: &[Message]) {
    enc.member("messages", 2);
    enc.begin_list(false);
    for (code, text) in messages {
        enc.element();
        enc.begin_sequence("Blast4-error");
        enc.member("code", 0);
        enc.int(*code);
        enc.end_member();
        enc.member("message", 1);
        enc.visible_str(text);
        enc.end_member();
        enc.end_sequence();
    }
    enc.end_list();
    enc.end_member();
}

fn write_ka_block<E: Encoder>(enc: &mut E, block: &KaBlock) {
    enc.element();
    enc.begin_sequence("Blast4-ka-block");
    enc.member("lambda", 0);
    enc.real(block.lambda);
    enc.end_member();
    enc.member("k", 1);
    enc.real(block.k);
    enc.end_member();
    enc.member("h", 2);
    enc.real(block.h);
    enc.end_member();
    enc.member("gapped", 3);
    enc.boolean(block.gapped);
    enc.end_member();
    enc.end_sequence();
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:352-364
/// Writes the `body` member of a `Blast4-request`: the `queue-search` request of one query
/// batch. The archive (`write_search_archive`) and the search strategy file
/// (`strategy::StrategyExport::write_request`) share it, as NCBI builds both with
/// `CExportStrategy`. `sets` is what the options handle holds beyond the command line's options
/// when the request is written (`OptionsHandle::request_sets`: the `EffectiveSearchSpace` option
/// the batches before left, RPS's database length and `.aux` costs, ...), applied onto the
/// request's options as NCBI's `x_AttachValue` applied them: a parameter of the command line is
/// replaced in its place (`-searchsp` keeps its place with the batch's value, `-searchsp 0`
/// too), the others follow in their order.
fn write_request_body<E: Encoder>(
    out: &mut ArchiveOut<'_, E>,
    request: &super::ArchiveRequest,
    batch: &RequestInputs<'_>,
    sets: &[Param],
    database_params: &[Param],
) {
    out.enc.member("body", 1);
    out.enc.variant("Blast4-request-body", "queue-search", 8);
    out.enc.begin_sequence("Blast4-queue-search-request");
    out.enc.member("program", 0);
    out.enc.visible_str(request.program.as_bytes());
    out.enc.end_member();
    out.enc.member("service", 1);
    out.enc.visible_str(request.service.as_bytes());
    out.enc.end_member();

    out.enc.member("queries", 2);
    match batch.queries {
        RequestQueries::Bioseqs(queries) => write_bioseq_queries(out, queries),
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:555-560
        // The object that the application read (`x_Process_Pssm` takes the application's `CRef`,
        // serialized at once), written with the request's encoder.
        RequestQueries::Pssm(pssm) => {
            out.enc.variant("Blast4-queries", "pssm", 0);
            crate::report::serial::pssm::write_value(&mut out.enc, pssm);
            out.enc.end_variant();
            out.flush_if_large();
        }
    }
    out.enc.end_member();

    out.enc.member("subject", 3);
    match &request.database {
        Some(name) => {
            out.enc.variant("Blast4-subject", "database", 0);
            out.enc.visible_str(name);
        }
        None => {
            out.enc.variant("Blast4-subject", "sequences", 1);
            out.enc.begin_list(false);
            for subject in batch.subjects {
                out.enc.element();
                write_bioseq(&mut out.enc, subject);
                out.flush_if_large();
            }
            out.enc.end_list();
        }
    }
    out.enc.end_variant();
    out.enc.end_member();

    write_request_options(out, request, batch, sets, database_params);
    out.enc.end_sequence();
    out.enc.end_variant();
    out.enc.end_member();
}

/// The `queries` of a request of FASTA queries: the batch's Seq-locs, or its Bioseqs.
fn write_bioseq_queries<E: Encoder>(out: &mut ArchiveOut<'_, E>, queries: &[BioseqInput<'_>]) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:609-617,648-655
    if queries.iter().all(|query| query.seq_loc.is_some()) {
        out.enc.variant("Blast4-queries", "seq-loc-list", 1);
        out.enc.begin_list(false);
        for loc in queries.iter().filter_map(|query| query.seq_loc.as_ref()) {
            out.enc.element();
            write_seq_loc(&mut out.enc, loc);
            out.flush_if_large();
        }
        out.enc.end_list();
        out.enc.end_variant();
    } else {
        out.enc.variant("Blast4-queries", "bioseq-set", 2);
        out.enc.begin_sequence("Bioseq-set");
        out.enc.member("seq-set", 7);
        out.enc.begin_list(false);
        for query in queries {
            out.enc.element();
            out.enc.variant("Seq-entry", "seq", 0);
            write_bioseq(&mut out.enc, query);
            out.enc.end_variant();
            out.flush_if_large();
        }
        out.enc.end_list();
        out.enc.end_member();
        out.enc.end_sequence();
        out.enc.end_variant();
    }
}

/// The options of a request: `algorithm-options` (the command line's with what the handle holds,
/// `sets`), `program-options` (the query's range and masks, the database's list and masking,
/// `database_params`) and `format-options` (psiblast's `Web_StepNumber`).
fn write_request_options<E: Encoder>(
    out: &mut ArchiveOut<'_, E>,
    request: &super::ArchiveRequest,
    batch: &RequestInputs<'_>,
    sets: &[Param],
    database_params: &[Param],
) {
    out.enc.member("algorithm-options", 5);
    out.enc.begin_named("Blast4-parameters");
    out.enc.begin_list(false);
    let set = |name: &str| sets.iter().find(|param| param.name() == name);
    for param in &request.algorithm {
        write_param(&mut out.enc, set(param.name()).unwrap_or(param));
    }
    for param in sets.iter().filter(|param| {
        !request
            .algorithm
            .iter()
            .any(|recorded| recorded.name() == param.name())
    }) {
        write_param(&mut out.enc, param);
    }
    out.enc.end_list();
    out.enc.end_named();
    out.enc.end_member();

    if batch.required.is_some() || !batch.lcase_masks.is_empty() || !database_params.is_empty() {
        out.enc.member("program-options", 6);
        out.enc.begin_named("Blast4-parameters");
        out.enc.begin_list(false);
        if let Some((start, end)) = batch.required {
            for (field, value) in [(fields::REQUIRED_START, start), (fields::REQUIRED_END, end)] {
                write_param(&mut out.enc, &Param::new(field, ParamValue::Integer(value)));
            }
        }
        for mask in batch.lcase_masks {
            out.enc.element();
            out.enc.begin_sequence("Blast4-parameter");
            out.enc.member("name", 0);
            out.enc.visible_str(fields::L_CASE_MASK.name().as_bytes());
            out.enc.end_member();
            out.enc.member("value", 1);
            let value_type = fields::L_CASE_MASK.value_type();
            out.enc
                .variant("Blast4-value", value_type.variant(), value_type.tag());
            write_mask(&mut out.enc, mask);
            out.enc.end_variant();
            out.enc.end_member();
            out.enc.end_sequence();
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:318-320
        // The database's list and masking come after the parameters of the query
        // (search_strategy.cpp:480-516: `TaxidList` or `NegativeTaxidList`, then
        // `DbFilteringAlgorithmKey` or `DbFilteringAlgorithmId`, then `SubjectMaskingType`).
        for param in database_params {
            write_param(&mut out.enc, param);
        }
        out.enc.end_list();
        out.enc.end_named();
        out.enc.end_member();
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:745-755
    if let Some(step) = request.step_number {
        out.enc.member("format-options", 7);
        out.enc.begin_named("Blast4-parameters");
        out.enc.begin_list(false);
        write_param(
            &mut out.enc,
            &Param::new(fields::WEB_STEP_NUMBER, ParamValue::Integer(step)),
        );
        out.enc.end_list();
        out.enc.end_named();
        out.enc.end_member();
    }
}

fn write_search_archive<E: Encoder, A: AlignSource>(
    out: &mut ArchiveOut<'_, E>,
    request: &super::ArchiveRequest,
    batch: &BatchArchive<'_>,
    sets: &[Param],
    database_params: &[Param],
    messages: &[Message],
    aligns: &mut A,
) -> Result<()> {
    out.enc.begin_document("Blast4-archive");
    out.enc.begin_sequence("Blast4-archive");
    out.enc.member("request", 0);
    out.enc.begin_sequence("Blast4-request");
    out.enc.member("ident", 0);
    // The request's `ident`: NCBI's `CBlastVersion::Print()`, LOSAT's version (`wording::LOSAT_VERSION`).
    out.enc
        .visible_str(crate::wording::LOSAT_VERSION.as_bytes());
    out.enc.end_member();
    write_request_body(out, request, &batch.request, sets, database_params);
    out.enc.end_sequence();
    out.enc.end_member();

    out.enc.member("results", 1);
    out.enc.begin_sequence("Blast4-get-search-results-reply");
    out.enc.member("alignments", 0);
    out.enc.begin_named("Seq-align-set");
    out.enc.begin_list(true);
    aligns.write(out)?;
    out.enc.end_list();
    out.enc.end_named();
    out.enc.end_member();

    out.enc.member("masks", 2);
    out.enc.begin_list(false);
    for mask in batch.masks {
        out.enc.element();
        write_mask(&mut out.enc, mask);
    }
    out.enc.end_list();
    out.enc.end_member();

    out.enc.member("ka-blocks", 3);
    out.enc.begin_list(false);
    for block in &batch.ka_blocks {
        write_ka_block(&mut out.enc, block);
    }
    out.enc.end_list();
    out.enc.end_member();

    out.enc.member("search-stats", 4);
    out.enc.begin_list(false);
    for text in [
        format!("Effective search space: {}", batch.search_space),
        format!("Effective search space used: {}", batch.search_space),
        format!("Length adjustment: {}", batch.length_adjustment),
    ] {
        out.enc.element();
        out.enc.visible_str(text.as_bytes());
    }
    out.enc.end_list();
    out.enc.end_member();
    out.enc.end_sequence();
    out.enc.end_member();

    if !messages.is_empty() {
        write_messages(&mut out.enc, messages);
    }
    out.enc.end_sequence();
    out.enc.end_document();
    Ok(())
}

/// Writes the `Blast4-request` of one query batch as a document of its own: the search strategy
/// file's (`-export_search_strategy`). NCBI's applications give `CExportStrategy` no client id
/// (blast_app_util.cpp:335-337), so the request has no `ident`.
fn write_strategy_document<E: Encoder>(
    out: &mut ArchiveOut<'_, E>,
    request: &super::ArchiveRequest,
    batch: &RequestInputs<'_>,
    sets: &[Param],
    database_params: &[Param],
) {
    out.enc.begin_document("Blast4-request");
    out.enc.begin_sequence("Blast4-request");
    write_request_body(out, request, batch, sets, database_params);
    out.enc.end_sequence();
    out.enc.end_document();
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:499-505
/// Whether `encoding` writes every string of a request's document as it is: the command line's
/// (`ArchiveRequest::written_as_is`), the batch's Bioseqs (`RequestInputs::written_as_is`) and the
/// string options that the run's handle adds (`sets`, and `database_params`: the
/// `DbFilteringAlgorithmKey` of a `-db_soft_mask` name, which the search strategy keeps after the
/// formatter's warning; NCBI 2.17.0 `blastn -db_soft_mask café`: `value string "caf##"` and
/// `Error: (803.7) [blastn] Blast4-request.body.queue-search.program-options.E.value.string`,
/// probe RG K/k4b/agent7/probe4 `bn_mask_e_f6`). The archive's and the search strategy's one rule
/// for writing through `Traced`.
fn request_written_as_is(
    encoding: ArchiveEncoding,
    request: &super::ArchiveRequest,
    batch: &RequestInputs<'_>,
    sets: &[Param],
    database_params: &[Param],
) -> bool {
    request.written_as_is(encoding)
        && batch.written_as_is(encoding)
        && encoding.writes_params(sets)
        && encoding.writes_params(database_params)
}

// No NCBI counterpart: carries an encoded `Blast4-request` together with its `COStreamBuffer`
// `BufferPiece` writes and `803.7` posts from batch encoding to `ostream::StrategyStream`.
/// One encoded `Blast4-request` of `-export_search_strategy` (`encode_strategy_request`):
/// the ASN.1 text bytes, NCBI's `COStreamBuffer` pieces (`BufferPiece`), NCBI's `803.7` posts
/// paired with the bytes `COStreamBuffer` had written before each post, and any `NullElement`
/// exception raised before the document completed.
pub(super) struct EncodedStrategyRequest {
    pub(super) bytes: Vec<u8>,
    pub(super) pieces: Vec<crate::report::serial::encoder::BufferPiece>,
    pub(super) posts: Vec<(u64, Vec<u8>)>,
    pub(super) thrown: Option<crate::cli::NativeError>,
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:635-651
// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:366-369
/// Encodes the search strategy's `Blast4-request` of one query batch (`write_strategy_document`)
/// in ASN.1 text through `Traced`, recording NCBI's `COStreamBuffer` pieces (`BufferPiece`) and
/// `803.7` posts (`Traced::take_posted`) so `ostream::StrategyStream` can feed each piece into
/// NCBI's 8,191-byte file buffer at the batch's start.
pub(super) fn encode_strategy_request(
    request: &super::ArchiveRequest,
    batch: &RequestInputs<'_>,
    sets: &[Param],
    database_params: &[Param],
) -> Result<EncodedStrategyRequest> {
    let mut bytes = Vec::new();
    let mut buf = Vec::new();
    let (pieces, posts, thrown) = {
        let mut out = ArchiveOut::new(Traced::new(AsnTextEncoder::new(&mut buf), 0), &mut bytes);
        write_strategy_document(&mut out, request, batch, sets, database_params);
        let thrown = out.enc.take_thrown();
        if thrown.is_some() {
            out.enc.flush_to(out.sink, true)?;
        } else if let Some(error) = out.enc.error() {
            return Err(error.into());
        }
        let posts = out
            .enc
            .take_posted()
            .into_iter()
            .map(|posted| {
                (
                    posted.written,
                    posted.replaced.post(Some(request.application)),
                )
            })
            .collect();
        if let Some(thrown) = thrown {
            (
                Vec::new(),
                posts,
                Some(crate::cli::NativeError::null_element(
                    request.application,
                    &thrown,
                )),
            )
        } else {
            let pieces = out.enc.pieces().to_vec();
            out.finish()?;
            (pieces, posts, None)
        }
    };
    Ok(EncodedStrategyRequest {
        bytes,
        pieces,
        posts,
        thrown,
    })
}

fn write_error_document<E: Encoder>(out: &mut ArchiveOut<'_, E>, messages: &[Message]) {
    out.enc.begin_document("Blast4-archive");
    out.enc.begin_sequence("Blast4-archive");
    out.enc.member("request", 0);
    out.enc.begin_sequence("Blast4-request");
    out.enc.member("body", 1);
    out.enc
        .variant("Blast4-request-body", "get-request-info", 9);
    out.enc.begin_sequence("Blast4-get-request-info-request");
    out.enc.member("request-id", 0);
    out.enc.visible_str(b"Error");
    out.enc.end_member();
    out.enc.end_sequence();
    out.enc.end_variant();
    out.enc.end_member();
    out.enc.end_sequence();
    out.enc.end_member();
    out.enc.member("results", 1);
    out.enc.begin_sequence("Blast4-get-search-results-reply");
    out.enc.end_sequence();
    out.enc.end_member();
    if !messages.is_empty() {
        write_messages(&mut out.enc, messages);
    }
    out.enc.end_sequence();
    out.enc.end_document();
}

/// Writes the archive of one query batch to `sink` and takes the messages saved since the
/// previous archive. `ARCHIVE_FORMAT` values that NCBI does not write leave the sink
/// untouched. A document with a string that the encoding does not write as it is
/// (`request_written_as_is`) is written through `Traced`,
/// and NCBI's posts of the bytes that its stream replaced go to `diagnostics` (`post_replaced`)
/// as NCBI writes them: each after the bytes of the document that NCBI's buffer had written to
/// the output by then (`Traced::posting`; NCBI's `cerr` is tied to `cout`).
///
/// The archive is one object of NCBI's object stream on the report's stream
/// (`ostream::Site::Serial`, `CBlastFormat::PrintArchive`: the archive inserted with `MSerial_AsnText`),
/// then `~CObjectOStream`'s flush. LOSAT inserts it in large parts; where one fails at the device,
/// NCBI's pieces of the archive (a second encoding of it through `Traced`, `archive_pieces`)
/// decide NCBI's end.
pub fn write_archive<A: AlignSource>(
    handle: &ArchiveHandle,
    batch: &BatchArchive<'_>,
    aligns: &mut A,
    sink: &mut dyn crate::report::ostream::ReportWrite,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    use crate::report::ostream::Site;
    let messages = handle.take_messages();
    sink.site(Site::Serial { whole: true })?;
    // The object and `Close`'s flush: where a device write of either fails, NCBI's pieces of the
    // object decide NCBI's end (LOSAT writes the archive in pieces of its own).
    let written = write_archive_object(handle, batch, aligns, &messages, sink, diagnostics)
        .and_then(|()| sink.ncbi_flush().map_err(Into::into));
    let undecided = written
        .as_ref()
        .err()
        .and_then(|error| error.downcast_ref::<std::io::Error>())
        .is_some_and(crate::report::ostream::is_undecided_report_write_failure);
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:265-276
    // NCBI's exception leaves `WriteArchive` before the application's `ResetMessages` (a failed
    // write of the device, or the throw of its stream at a `NULL` list element,
    // `encoder::NullElement`): the archive's messages are saved still, in their order, for the
    // error archive (`PrintErrorArchive` of the application's `CATCH_ALL`).
    let failed = |error: anyhow::Error| {
        handle.restore_messages(messages.clone());
        error
    };
    if undecided {
        let pieces = archive_pieces(handle, batch, aligns, &messages)?;
        sink.resolve_object(&pieces)
            .map_err(|error| failed(error.into()))?;
    } else {
        written.map_err(failed)?;
    }
    sink.site(Site::Report)?;
    Ok(())
}

/// NCBI's insertions of the archive of `batch` into the report's stream (`Traced::pieces`): the
/// writes of its object stream's buffer, from a second encoding of the archive (with the
/// `messages` of the first), whose bytes are dropped.
#[cold]
#[inline(never)]
fn archive_pieces<A: AlignSource>(
    handle: &ArchiveHandle,
    batch: &BatchArchive<'_>,
    aligns: &mut A,
    messages: &[Message],
) -> Result<Vec<usize>> {
    let sets = handle.request_sets(batch.searched_space);
    let params = handle.options().database_params();
    let write = SearchArchive {
        request: handle.request(),
        batch,
        sets: &sets,
        database_params: &params,
        messages,
    };
    let mut buf = Vec::new();
    let mut dropped = std::io::sink();
    Ok(match handle.encoding() {
        ArchiveEncoding::Nothing => Vec::new(),
        ArchiveEncoding::Xml => write
            .to(
                Traced::xml(XmlEncoder::archive(&mut buf), 0),
                &mut dropped,
                aligns,
            )?
            .enc
            .piece_sizes(),
        ArchiveEncoding::Text => write
            .to(
                Traced::new(AsnTextEncoder::new(&mut buf), 0),
                &mut dropped,
                aligns,
            )?
            .enc
            .piece_sizes(),
        ArchiveEncoding::Binary => write
            .to(
                Traced::new(AsnBinaryEncoder::new(&mut buf), 0),
                &mut dropped,
                aligns,
            )?
            .enc
            .piece_sizes(),
    })
}

/// The archive of `write_archive`, written through the encoder of the handle's encoding.
fn write_archive_object<A: AlignSource>(
    handle: &ArchiveHandle,
    batch: &BatchArchive<'_>,
    aligns: &mut A,
    messages: &[Message],
    sink: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    let mut sets = handle.request_sets(batch.searched_space);
    let round_request;
    let request = match batch.round {
        Some(round) => {
            if let Some(mode) = round.comp_based_stats {
                sets.retain(|param| param.field() != fields::COMPOSITION_BASED_STATS);
                sets.push(Param::new(
                    fields::COMPOSITION_BASED_STATS,
                    ParamValue::Integer(mode),
                ));
            }
            let request = handle.request();
            round_request = super::ArchiveRequest {
                step_number: request.database.is_some().then_some(round.number),
                ..request.clone()
            };
            &round_request
        }
        None => handle.request(),
    };
    let encoding = handle.encoding();
    // A big integer of -2^63 (`-dbsize`), which NCBI's ASN.1 text and XML writers print as `-0` (a
    // signed overflow of `COStreamBuffer::PutInt8`), is written as its digits: PD-LOSAT-NCBI-DEFECTS
    // approved exception 17, Int8 -2^63 written as digits (the Owner's Q1 = A of 2026-11-17,
    // `big_int`).
    let params = handle.options().database_params();
    let as_is = request_written_as_is(encoding, request, &batch.request, &sets, &params);
    let post = |replaced: &ReplacedChar| {
        super::post_replaced(
            diagnostics,
            handle,
            request.application,
            std::slice::from_ref(replaced),
        )?;
        diagnostics.flush()
    };
    let mut buf = Vec::new();
    let write = SearchArchive {
        request,
        batch,
        sets: &sets,
        database_params: &params,
        messages,
    };
    let application = request.application;
    let written = match (encoding, as_is) {
        (ArchiveEncoding::Nothing, _) => Ok(()),
        (ArchiveEncoding::Xml, true) => write
            .to(XmlEncoder::archive(&mut buf), sink, aligns)?
            .finish(),
        (ArchiveEncoding::Xml, false) => write
            .to(
                Traced::xml(XmlEncoder::archive(&mut buf), 0).posting(post),
                sink,
                aligns,
            )?
            .finish_archive(application),
        (ArchiveEncoding::Text, true) => write
            .to(AsnTextEncoder::new(&mut buf), sink, aligns)?
            .finish(),
        (ArchiveEncoding::Text, false) => write
            .to(
                Traced::new(AsnTextEncoder::new(&mut buf), 0).posting(post),
                sink,
                aligns,
            )?
            .finish_archive(application),
        (ArchiveEncoding::Binary, true) => write
            .to(AsnBinaryEncoder::new(&mut buf), sink, aligns)?
            .finish(),
        (ArchiveEncoding::Binary, false) => write
            .to(
                Traced::new(AsnBinaryEncoder::new(&mut buf), 0).posting(post),
                sink,
                aligns,
            )?
            .finish_archive(application),
    };
    written
}

/// What `write_search_archive` writes besides the encoder and the Seq-aligns.
struct SearchArchive<'a> {
    request: &'a super::ArchiveRequest,
    batch: &'a BatchArchive<'a>,
    sets: &'a [Param],
    database_params: &'a [Param],
    messages: &'a [Message],
}

impl SearchArchive<'_> {
    /// The document written through `enc` to `sink`, not finished.
    fn to<'s, E: Encoder, A: AlignSource>(
        &self,
        enc: E,
        sink: &'s mut dyn Write,
        aligns: &mut A,
    ) -> Result<ArchiveOut<'s, E>> {
        let mut out = ArchiveOut::new(enc, sink);
        write_search_archive(
            &mut out,
            self.request,
            self.batch,
            self.sets,
            self.database_params,
            self.messages,
            aligns,
        )?;
        Ok(out)
    }
}

/// Writes the error archive (`PrintErrorArchive`): the messages saved since the last
/// archive and `messages` (the error lines). Nothing is written when the command line's
/// `-outfmt` text is not the integer 11 (`writes_error_archive`).
///
/// The archive is one object on the report's stream (`ostream::Site::Serial`); written to a
/// stream that failed in the run (NCBI's `PrintErrorArchive` after `CATCH_ALL`, the exception mask
/// reset), its buffer's first write fails: NCBI's write fault, posted at that write's position
/// (`error_archive_fault`), and nothing is written.
pub fn write_error_archive(
    handle: &ArchiveHandle,
    error_lines: &[Message],
    sink: &mut dyn crate::report::ostream::ReportWrite,
) -> Result<()> {
    use crate::report::ostream::Site;
    if !handle.error_archive_enabled() {
        return Ok(());
    }
    let mut messages = handle.take_messages();
    messages.extend(error_lines.iter().cloned());
    sink.site(Site::Serial { whole: true })?;
    let written = write_error_object(handle, &messages, sink);
    sink.site(Site::Report)?;
    match written {
        Err(error)
            if error
                .downcast_ref::<std::io::Error>()
                .is_some_and(crate::report::ostream::is_undecided_report_write_failure) =>
        {
            Ok(sink.write_fault(error_archive_fault(handle.encoding(), &messages))?)
        }
        written => written,
    }
}

/// Where NCBI's buffer of the error archive of `messages` first writes to its output (the position
/// of its write fault: `OutBuffer::first_write`, the bytes it had put into the buffer; for a small
/// archive the document without the ASN.1 text's separator, which is a write of its own).
#[cold]
#[inline(never)]
pub(super) fn error_archive_fault(encoding: ArchiveEncoding, messages: &[Message]) -> usize {
    let mut buf = Vec::new();
    let mut dropped = std::io::sink();
    let first = match encoding {
        ArchiveEncoding::Nothing => return 0,
        ArchiveEncoding::Xml => {
            let mut out =
                ArchiveOut::new(Traced::xml(XmlEncoder::archive(&mut buf), 0), &mut dropped);
            write_error_document(&mut out, messages);
            out.enc.first_write()
        }
        ArchiveEncoding::Text => {
            let mut out =
                ArchiveOut::new(Traced::new(AsnTextEncoder::new(&mut buf), 0), &mut dropped);
            write_error_document(&mut out, messages);
            out.enc.first_write()
        }
        ArchiveEncoding::Binary => {
            let mut out = ArchiveOut::new(
                Traced::new(AsnBinaryEncoder::new(&mut buf), 0),
                &mut dropped,
            );
            write_error_document(&mut out, messages);
            out.enc.first_write()
        }
    };
    first.unwrap_or(0) as usize
}

/// The error archive of `messages` written through the encoder of the handle's encoding.
fn write_error_object(
    handle: &ArchiveHandle,
    messages: &[Message],
    sink: &mut dyn Write,
) -> Result<()> {
    match handle.encoding() {
        ArchiveEncoding::Nothing => Ok(()),
        ArchiveEncoding::Xml => {
            let mut buf = Vec::new();
            let mut out = ArchiveOut::new(XmlEncoder::archive(&mut buf), sink);
            write_error_document(&mut out, messages);
            out.finish()
        }
        ArchiveEncoding::Text => {
            let mut buf = Vec::new();
            let mut out = ArchiveOut::new(AsnTextEncoder::new(&mut buf), sink);
            write_error_document(&mut out, messages);
            out.finish()
        }
        ArchiveEncoding::Binary => {
            let mut buf = Vec::new();
            let mut out = ArchiveOut::new(AsnBinaryEncoder::new(&mut buf), sink);
            write_error_document(&mut out, messages);
            out.finish()
        }
    }
}
