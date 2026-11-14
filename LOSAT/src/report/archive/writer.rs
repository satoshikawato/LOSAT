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
use super::masks::{write_mask, MaskEntry};
use super::options::{Param, ParamValue};
use super::{ArchiveEncoding, ArchiveHandle};
use crate::report::diag_log::Message;
use crate::report::serial::asn_ber::AsnBinaryEncoder;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::encoder::{Encoder, ReplacedChar, Traced};
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

/// What the request of one query batch reads besides the run's options (`ArchiveRequest`): the
/// archive's and the search strategy's (`requests::RunRequests` makes it from the batch's
/// records).
#[derive(Clone, Copy)]
pub struct RequestInputs<'a> {
    /// The batch's queries, in order.
    pub queries: &'a [BioseqInput<'a>],
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
    /// its reader refuses or replaces where NCBI would post).
    pub fn written_as_is(&self, encoding: ArchiveEncoding) -> bool {
        self.queries
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
    /// The batch's `EffectiveSearchSpace` option when the batch was searched.
    pub searched_space: Option<i64>,
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

impl<E: Encoder> ArchiveOut<'_, Traced<E>> {
    /// `finish` of a document written through `Traced` without a post function (the search
    /// strategy's, whose file is apart from the output), after `post` of all the bytes that its
    /// stream replaced (what was handed to the sink goes first). Another error of the encoder
    /// posts nothing.
    fn finish_posting(
        mut self,
        post: impl FnOnce(&[ReplacedChar]) -> std::io::Result<()>,
    ) -> Result<()> {
        if self.enc.error().is_none() {
            let replaced = self.enc.take_replaced_chars();
            if !replaced.is_empty() {
                if self.io_error.is_none() {
                    if let Err(error) = self.sink.flush() {
                        self.io_error = Some(error);
                    }
                }
                post(&replaced)?;
            }
        }
        self.finish()
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
    enc.visible_str(param.name.as_bytes());
    enc.end_member();
    enc.member("value", 1);
    match &param.value {
        ParamValue::Evalue(value) => {
            enc.variant("Blast4-value", "cutoff", 3);
            enc.variant("Blast4-cutoff", "e-value", 0);
            enc.real(*value);
            enc.end_variant();
        }
        ParamValue::Integer(value) => {
            enc.variant("Blast4-value", "integer", 4);
            enc.int(*value);
        }
        ParamValue::BigInteger(value) => {
            enc.variant("Blast4-value", "big-integer", 0);
            enc.big_int(*value);
        }
        ParamValue::Boolean(value) => {
            enc.variant("Blast4-value", "boolean", 2);
            enc.boolean(*value);
        }
        ParamValue::Real(value) => {
            enc.variant("Blast4-value", "real", 6);
            enc.real(*value);
        }
        ParamValue::Str(value) => {
            enc.variant("Blast4-value", "string", 11);
            enc.visible_str(value.as_bytes());
        }
        ParamValue::IntegerList(values) => {
            enc.variant("Blast4-value", "integer-list", 16);
            enc.begin_list(false);
            for &value in values {
                enc.element();
                enc.int(value);
            }
            enc.end_list();
        }
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
// ```c++
// CRef<objects::CBlast4_request> CExportStrategy::GetSearchStrategy(void)
// {
// 	CRef<CBlast4_request> retval(new CBlast4_request);
// 	if (!m_ClientId.empty())
// 	{
//         retval->SetIdent(m_ClientId);
//     }
//
//     CRef<CBlast4_request_body> body(new CBlast4_request_body);
//     body->SetQueue_search(*m_QueueSearchRequest);
//     retval->SetBody(*body);
//     return retval;
// }
// ```
/// Writes the `body` member of a `Blast4-request`: the `queue-search` request of one query
/// batch. The archive (`write_search_archive`) and the search strategy file
/// (`strategy::StrategyExport::write_request`) share it, as NCBI builds both with
/// `CExportStrategy`. `search_space_option` is the `EffectiveSearchSpace` option that the
/// request has when it is written.
fn write_request_body<E: Encoder>(
    out: &mut ArchiveOut<'_, E>,
    request: &super::ArchiveRequest,
    batch: &RequestInputs<'_>,
    search_space_option: Option<i64>,
    database_list: Option<&Param>,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:609-617,648-655
    // ```c++
    //         ITERATE(IRemoteQueryData::TSeqLocs, itr, seqloc_list)
    //         {
    //             if (IsLocalId((*itr)->GetId()))
    //             {
    //                 has_local_ids = true;
    //                 break;
    //             }
    //         }
    //     ...
    //     if (has_local_ids)
    //     {
    //         Q->SetBioseq_set(*bioseq_set);
    //     }
    //     else
    //     {
    //     	Q->SetSeq_loc_list() = seqloc_list;
    //     }
    // ```
    if batch.queries.iter().all(|query| query.seq_loc.is_some()) {
        out.enc.variant("Blast4-queries", "seq-loc-list", 1);
        out.enc.begin_list(false);
        for loc in batch
            .queries
            .iter()
            .filter_map(|query| query.seq_loc.as_ref())
        {
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
        for query in batch.queries {
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

    out.enc.member("algorithm-options", 5);
    out.enc.begin_named("Blast4-parameters");
    out.enc.begin_list(false);
    let mut written_search_space = false;
    for param in &request.algorithm {
        // `-searchsp` holds the option's place; its value is the batch's (`-searchsp 0` too).
        if let ("EffectiveSearchSpace", Some(space)) = (param.name, search_space_option) {
            written_search_space = true;
            write_param(
                &mut out.enc,
                &Param {
                    name: param.name,
                    value: ParamValue::BigInteger(space),
                },
            );
            continue;
        }
        written_search_space |= param.name == "EffectiveSearchSpace";
        write_param(&mut out.enc, param);
    }
    if let (Some(space), false) = (search_space_option, written_search_space) {
        write_param(
            &mut out.enc,
            &Param {
                name: "EffectiveSearchSpace",
                value: ParamValue::BigInteger(space),
            },
        );
    }
    out.enc.end_list();
    out.enc.end_named();
    out.enc.end_member();

    if batch.required.is_some() || !batch.lcase_masks.is_empty() || database_list.is_some() {
        out.enc.member("program-options", 6);
        out.enc.begin_named("Blast4-parameters");
        out.enc.begin_list(false);
        if let Some((start, end)) = batch.required {
            for (name, value) in [("RequiredStart", start), ("RequiredEnd", end)] {
                write_param(
                    &mut out.enc,
                    &Param {
                        name,
                        value: ParamValue::Integer(value),
                    },
                );
            }
        }
        for mask in batch.lcase_masks {
            out.enc.element();
            out.enc.begin_sequence("Blast4-parameter");
            out.enc.member("name", 0);
            out.enc.visible_str(b"LCaseMask");
            out.enc.end_member();
            out.enc.member("value", 1);
            out.enc.variant("Blast4-value", "query-mask", 26);
            write_mask(&mut out.enc, mask);
            out.enc.end_variant();
            out.enc.end_member();
            out.enc.end_sequence();
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:318-320
        // ```c++
        //     x_Process_BlastOptions(opts_handle);
        //     x_Process_Query(query);
        //     x_Process_SearchDb(db);
        // ```
        // The database's list comes after the parameters of the query.
        if let Some(param) = database_list {
            write_param(&mut out.enc, param);
        }
        out.enc.end_list();
        out.enc.end_named();
        out.enc.end_member();
    }
    out.enc.end_sequence();
    out.enc.end_variant();
    out.enc.end_member();
}

fn write_search_archive<E: Encoder, A: AlignSource>(
    out: &mut ArchiveOut<'_, E>,
    request: &super::ArchiveRequest,
    batch: &BatchArchive<'_>,
    search_space_option: Option<i64>,
    database_list: Option<&Param>,
    messages: &[Message],
    aligns: &mut A,
) -> Result<()> {
    out.enc.begin_document("Blast4-archive");
    out.enc.begin_sequence("Blast4-archive");
    out.enc.member("request", 0);
    out.enc.begin_sequence("Blast4-request");
    out.enc.member("ident", 0);
    // The request's `ident` is `CBlastVersion::Print()` (`crate::NCBI_BLAST_VERSION`).
    out.enc.visible_str(crate::NCBI_BLAST_VERSION.as_bytes());
    out.enc.end_member();
    write_request_body(
        out,
        request,
        &batch.request,
        search_space_option,
        database_list,
    );
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
    search_space_option: Option<i64>,
    database_list: Option<&Param>,
) {
    out.enc.begin_document("Blast4-request");
    out.enc.begin_sequence("Blast4-request");
    write_request_body(out, request, batch, search_space_option, database_list);
    out.enc.end_sequence();
    out.enc.end_document();
}

/// Writes the search strategy's request of one query batch (`write_strategy_document`) to `sink`
/// in ASN.1 text (the strategy file is always ASN.1 text), through `Traced` when ASN.1 text does
/// not write one of its strings as it is (`ArchiveRequest::written_as_is`,
/// `RequestInputs::written_as_is`): `post` then gets the bytes that it replaced.
pub(super) fn write_strategy_request(
    sink: &mut dyn Write,
    request: &super::ArchiveRequest,
    batch: &RequestInputs<'_>,
    search_space_option: Option<i64>,
    database_list: Option<&Param>,
    post: impl FnOnce(&[ReplacedChar]) -> std::io::Result<()>,
) -> Result<()> {
    let mut buf = Vec::new();
    let encoding = ArchiveEncoding::Text;
    if request.written_as_is(encoding) && batch.written_as_is(encoding) {
        let mut out = ArchiveOut::new(AsnTextEncoder::new(&mut buf), sink);
        write_strategy_document(&mut out, request, batch, search_space_option, database_list);
        out.finish()
    } else {
        let mut out = ArchiveOut::new(Traced::new(AsnTextEncoder::new(&mut buf), 0), sink);
        write_strategy_document(&mut out, request, batch, search_space_option, database_list);
        out.finish_posting(post)
    }
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
/// (`ArchiveRequest::written_as_is`, `RequestInputs::written_as_is`) is written through `Traced`,
/// and NCBI's posts of the bytes that its stream replaced go to `diagnostics` (`post_replaced`)
/// as NCBI writes them: each after the bytes of the document that NCBI's buffer had written to
/// the output by then (`Traced::posting`; NCBI's `cerr` is tied to `cout`).
pub fn write_archive<A: AlignSource>(
    handle: &ArchiveHandle,
    batch: &BatchArchive<'_>,
    aligns: &mut A,
    sink: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    let messages = handle.take_messages();
    let space = handle.search_space_option(batch.searched_space);
    let list = handle.options().database_list();
    let request = handle.request();
    let encoding = handle.encoding();
    let as_is = request.written_as_is(encoding) && batch.request.written_as_is(encoding);
    let post = |replaced: &ReplacedChar| {
        super::post_replaced(
            diagnostics,
            Some(handle),
            false,
            request.program,
            std::slice::from_ref(replaced),
        )?;
        diagnostics.flush()
    };
    let mut buf = Vec::new();
    let write = SearchArchive {
        request,
        batch,
        search_space_option: space,
        database_list: list,
        messages: &messages,
    };
    match (encoding, as_is) {
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
            .finish(),
        (ArchiveEncoding::Text, true) => write
            .to(AsnTextEncoder::new(&mut buf), sink, aligns)?
            .finish(),
        (ArchiveEncoding::Text, false) => write
            .to(
                Traced::new(AsnTextEncoder::new(&mut buf), 0).posting(post),
                sink,
                aligns,
            )?
            .finish(),
        (ArchiveEncoding::Binary, true) => write
            .to(AsnBinaryEncoder::new(&mut buf), sink, aligns)?
            .finish(),
        (ArchiveEncoding::Binary, false) => write
            .to(
                Traced::new(AsnBinaryEncoder::new(&mut buf), 0).posting(post),
                sink,
                aligns,
            )?
            .finish(),
    }
}

/// What `write_search_archive` writes besides the encoder and the Seq-aligns.
struct SearchArchive<'a> {
    request: &'a super::ArchiveRequest,
    batch: &'a BatchArchive<'a>,
    search_space_option: Option<i64>,
    database_list: Option<&'a Param>,
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
            self.search_space_option,
            self.database_list,
            self.messages,
            aligns,
        )?;
        Ok(out)
    }
}

/// Writes the error archive (`PrintErrorArchive`): the messages saved since the last
/// archive and `messages` (the error lines). Nothing is written when the command line's
/// `-outfmt` text is not the integer 11 (`writes_error_archive`).
pub fn write_error_archive(
    handle: &ArchiveHandle,
    error_lines: &[Message],
    sink: &mut dyn Write,
) -> Result<()> {
    if !handle.error_archive_enabled() {
        return Ok(());
    }
    let mut messages = handle.take_messages();
    messages.extend(error_lines.iter().cloned());
    match handle.encoding() {
        ArchiveEncoding::Nothing => Ok(()),
        ArchiveEncoding::Xml => {
            let mut buf = Vec::new();
            let mut out = ArchiveOut::new(XmlEncoder::archive(&mut buf), sink);
            write_error_document(&mut out, &messages);
            out.finish()
        }
        ArchiveEncoding::Text => {
            let mut buf = Vec::new();
            let mut out = ArchiveOut::new(AsnTextEncoder::new(&mut buf), sink);
            write_error_document(&mut out, &messages);
            out.finish()
        }
        ArchiveEncoding::Binary => {
            let mut buf = Vec::new();
            let mut out = ArchiveOut::new(AsnBinaryEncoder::new(&mut buf), sink);
            write_error_document(&mut out, &messages);
            out.finish()
        }
    }
}
