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
use crate::report::serial::encoder::{Encoder, SerialError, Traced};
use crate::report::serial::visible::BadChars;
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
/// do not return one. A document is written with a plain encoder first (`ArchiveOut::plain`), which
/// hands to the sink only the bytes older than NCBI's buffer can hold (`HOLD`); when the encoder
/// reports a byte that NCBI writes as `#` (`SerialError::NonVisibleChar`), it hands out nothing
/// more and `finish` asks for the document again (`Finish::Again`), which is then written through
/// `Traced` and NCBI's form of the encoder (`ArchiveOut::traced`, `write_asn`), from the byte where
/// the first writing stopped handing out. With `posts`, the messages of the bytes the document
/// writes as `#` go to the run's diagnostics among the document's bytes, where NCBI's buffer had
/// handed them out (`flush`).
pub struct ArchiveOut<'s, E: Encoder> {
    pub enc: E,
    sink: &'s mut dyn Write,
    io_error: Option<std::io::Error>,
    posts: Option<PostOut<'s>>,
    /// The bytes that a flush keeps (`HOLD` in the plain writing, 0 through `Traced`).
    keep: usize,
    /// The plain writing: the bytes handed to the sink.
    handed: usize,
    /// Through `Traced` without `posts`: the strings with replaced bytes, for the caller.
    bad_chars: Vec<BadChars>,
    /// Through `Traced`: the bytes of the document still to drop instead of handing them out
    /// (those that the first writing handed out).
    skip: usize,
}

/// The bytes that the plain writing of a document keeps from the sink: NCBI's buffer never holds
/// more (`stream_buffer::INITIAL_SIZE`: none of the archive's reservations needs a larger buffer,
/// and `Write` fills it without growing it), so the message of a byte found later falls after
/// every byte handed out (NCBI has handed out at least the document's bytes but its buffer's).
const HOLD: usize = crate::report::serial::stream_buffer::INITIAL_SIZE;

/// How the writing of a document ended (`ArchiveOut::finish`).
pub(super) enum Finish {
    /// The document is written, with the strings that the encoder wrote with bytes replaced by
    /// `#` and has not posted (`ArchiveOut::traced` without `posts`).
    Done(Vec<BadChars>),
    /// The plain encoder reported a byte that NCBI writes as `#`: write the document again
    /// through `Traced`, from the byte (the bytes handed out).
    Again(usize),
}

/// Where an archive's messages go (`ArchiveHandle::write_unsaved_text`).
pub(super) struct PostOut<'s> {
    pub(super) diagnostics: &'s mut dyn Write,
    pub(super) handle: &'s ArchiveHandle,
}

/// Hands the encoder's final bytes to `sink` (`all`: every byte), with the messages of the bytes
/// written as `#` so far, each after the document's bytes that NCBI's buffer has handed to its
/// output when it posts it (`visible::posts_after`: on a stream that both go to, NCBI's message
/// falls there). A plain encoder that reported such a byte hands nothing more out: its bytes are
/// dropped (the document is written again, `Finish::Again`).
#[allow(clippy::too_many_arguments)]
fn flush<E: Encoder>(
    enc: &mut E,
    sink: &mut dyn Write,
    posts: &mut Option<PostOut<'_>>,
    bad_chars: &mut Vec<BadChars>,
    keep: usize,
    handed: &mut usize,
    skip: &mut usize,
    all: bool,
) -> std::io::Result<()> {
    if *skip > 0 {
        let mut skipping = SkipFirst {
            inner: sink,
            left: *skip,
        };
        let flushed = flush(
            enc,
            &mut skipping,
            posts,
            bad_chars,
            keep,
            handed,
            &mut 0,
            all,
        );
        *skip = skipping.left;
        return flushed;
    }
    if keep > 0 {
        if let Some(SerialError::NonVisibleChar(_)) = enc.error() {
            return enc.flush_keeping(&mut std::io::sink(), keep);
        }
        let before = enc.buffered();
        if all {
            enc.flush_to(sink, true)?;
        } else {
            enc.flush_keeping(sink, keep)?;
        }
        *handed += before - enc.buffered();
        return Ok(());
    }
    match posts {
        Some(out) => {
            let program = Some(out.handle.request().program);
            for (after, text) in
                crate::report::serial::visible::posts_after(&enc.take_bad_chars(), program)
            {
                enc.flush_through(sink, after)?;
                sink.flush()?;
                out.handle
                    .write_unsaved_text(&mut *out.diagnostics, &text)?;
            }
        }
        None => bad_chars.extend(enc.take_bad_chars()),
    }
    enc.flush_to(sink, all)
}

impl<'s, E: Encoder> ArchiveOut<'s, E> {
    /// The first writing of a document, with a plain encoder.
    pub(super) fn plain(enc: E, sink: &'s mut dyn Write) -> Self {
        Self::new(enc, sink, None, HOLD)
    }

    /// A document written once, whose encoder writes no `#` (the XML archive, the error archive
    /// of sanitized messages).
    pub(super) fn once(enc: E, sink: &'s mut dyn Write) -> Self {
        Self::new(enc, sink, None, 0)
    }

    /// The second writing of a document, through `Traced` and NCBI's form of an encoder; the
    /// messages go to `posts` among the bytes, or to the caller (`Finish::Done`).
    /// (`skip`: the bytes of the document that the first writing handed out.)
    pub(super) fn traced(
        enc: E,
        sink: &'s mut dyn Write,
        posts: Option<PostOut<'s>>,
        skip: usize,
    ) -> Self {
        let mut out = Self::new(enc, sink, posts, 0);
        out.skip = skip;
        out
    }

    fn new(enc: E, sink: &'s mut dyn Write, posts: Option<PostOut<'s>>, keep: usize) -> Self {
        ArchiveOut {
            enc,
            sink,
            io_error: None,
            posts,
            keep,
            handed: 0,
            bad_chars: Vec::new(),
            skip: 0,
        }
    }

    /// Hands the encoder's final bytes to the sink once enough have gathered. Call it after
    /// each element of a long list.
    pub fn flush_if_large(&mut self) {
        if self.io_error.is_none() && self.enc.buffered() >= FLUSH_BYTES + self.keep {
            if let Err(error) = flush(
                &mut self.enc,
                self.sink,
                &mut self.posts,
                &mut self.bad_chars,
                self.keep,
                &mut self.handed,
                &mut self.skip,
                false,
            ) {
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
            posts,
            keep,
            handed,
            bad_chars,
            skip,
        } = self;
        let keep = *keep;
        let mut flush_large = |enc: &mut E| {
            if io_error.is_none() && enc.buffered() >= FLUSH_BYTES + keep {
                if let Err(error) = flush(
                    enc,
                    &mut **sink,
                    posts,
                    bad_chars,
                    keep,
                    handed,
                    skip,
                    false,
                ) {
                    *io_error = Some(error);
                }
            }
        };
        write(enc, &mut flush_large)
    }

    /// Hands the rest of the document to the sink, or asks for it again (`Finish::Again`) when
    /// the plain encoder reported a byte that NCBI writes as `#`.
    pub(super) fn finish(mut self) -> Result<Finish> {
        if self.keep > 0 {
            if let Some(SerialError::NonVisibleChar(_)) = self.enc.error() {
                if let Some(error) = self.io_error {
                    return Err(error.into());
                }
                return Ok(Finish::Again(self.handed));
            }
        }
        if self.io_error.is_none() {
            if let Err(error) = flush(
                &mut self.enc,
                self.sink,
                &mut self.posts,
                &mut self.bad_chars,
                self.keep,
                &mut self.handed,
                &mut self.skip,
                true,
            ) {
                self.io_error = Some(error);
            }
        }
        if let Some(error) = self.enc.error() {
            return Err(error.into());
        }
        match self.io_error {
            Some(error) => Err(error.into()),
            None => Ok(Finish::Done(self.bad_chars)),
        }
    }
}

/// A document of the archive's writer, written with any encoder (`write_asn`).
pub(super) trait ArchiveDocument {
    fn write<E: Encoder>(&mut self, out: &mut ArchiveOut<'_, E>) -> Result<()>;
}

/// Writes `document` to `sink` in ASN.1 text or binary (`binary`): once with the plain encoder,
/// and when it reported a byte that NCBI writes as `#`, again through `Traced` (`write_again`).
/// The strings with replaced bytes whose messages did not go to `posts`.
pub(super) fn write_asn<'s, D: ArchiveDocument>(
    binary: bool,
    document: &mut D,
    sink: &'s mut dyn Write,
    posts: Option<PostOut<'s>>,
) -> Result<Vec<BadChars>> {
    let mut buf = Vec::new();
    let finish = if binary {
        let mut out = ArchiveOut::plain(AsnBinaryEncoder::new(&mut buf), &mut *sink);
        document.write(&mut out)?;
        out.finish()?
    } else {
        let mut out = ArchiveOut::plain(AsnTextEncoder::new(&mut buf), &mut *sink);
        document.write(&mut out)?;
        out.finish()?
    };
    match finish {
        Finish::Done(bad_chars) => Ok(bad_chars),
        Finish::Again(handed) => write_again(binary, document, sink, posts, handed),
    }
}

/// `document` written again through `Traced` of NCBI's form of the encoder, its first `handed`
/// bytes (which the first writing handed out) not handed out again.
#[cold]
#[inline(never)]
fn write_again<'s, D: ArchiveDocument>(
    binary: bool,
    document: &mut D,
    sink: &'s mut dyn Write,
    posts: Option<PostOut<'s>>,
    handed: usize,
) -> Result<Vec<BadChars>> {
    let mut buf = Vec::new();
    let finish = if binary {
        let enc = Traced::new(AsnBinaryEncoder::ncbi(&mut buf));
        let mut out = ArchiveOut::traced(enc, sink, posts, handed);
        document.write(&mut out)?;
        out.finish()?
    } else {
        let enc = Traced::new(AsnTextEncoder::ncbi(&mut buf));
        let mut out = ArchiveOut::traced(enc, sink, posts, handed);
        document.write(&mut out)?;
        out.finish()?
    };
    match finish {
        Finish::Done(bad_chars) => Ok(bad_chars),
        Finish::Again(_) => unreachable!("only the plain writing asks for the document again"),
    }
}

/// A sink that drops the first `left` bytes written to it (the bytes of a document written again
/// that the first writing handed out).
struct SkipFirst<'w> {
    inner: &'w mut dyn Write,
    left: usize,
}

impl Write for SkipFirst<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let skipped = self.left.min(buf.len());
        self.left -= skipped;
        self.inner.write_all(&buf[skipped..])?;
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
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
pub(super) fn write_strategy_document<E: Encoder>(
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

/// The `Blast4-request` of one query batch (`write_strategy_document`) as an `ArchiveDocument`.
pub(super) struct StrategyDocument<'h> {
    pub(super) request: &'h super::ArchiveRequest,
    pub(super) batch: &'h RequestInputs<'h>,
    pub(super) search_space_option: Option<i64>,
    pub(super) database_list: Option<&'h Param>,
}

impl ArchiveDocument for StrategyDocument<'_> {
    fn write<E: Encoder>(&mut self, out: &mut ArchiveOut<'_, E>) -> Result<()> {
        write_strategy_document(
            out,
            self.request,
            self.batch,
            self.search_space_option,
            self.database_list,
        );
        Ok(())
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

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:304-306
// ```c++
//             if (isArchiveFormat) {
//                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
//                 m_Bah.ResetMessages();
// ```
/// Writes the archive of one query batch to `sink` and takes the messages saved since the
/// previous archive. `ARCHIVE_FORMAT` values that NCBI does not write leave the sink
/// untouched. The messages of the bytes that the ASN.1 text or binary stream writes as `#` (a
/// title, the `-db` name or the `-window_masker_db` path outside 0x20..0x7E) go to `diagnostics`
/// as NCBI posts them while it writes the archive, and are not saved for the next archive:
/// NCBI's handler saves them, and `ResetMessages` after `WriteArchive` drops them
/// (`ArchiveHandle::write_unsaved_posts`).
pub fn write_archive<A: AlignSource>(
    handle: &ArchiveHandle,
    batch: &BatchArchive<'_>,
    aligns: &mut A,
    sink: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    let messages = handle.take_messages_through(batch.request.queries.len());
    let mut document = SearchArchive {
        request: handle.request(),
        batch,
        search_space_option: handle.search_space_option(batch.searched_space),
        database_list: handle.options().database_list(),
        messages: &messages,
        aligns,
    };
    let binary = match handle.encoding() {
        ArchiveEncoding::Nothing => return Ok(()),
        ArchiveEncoding::Xml => {
            let mut buf = Vec::new();
            let mut out = ArchiveOut::once(XmlEncoder::archive(&mut buf), sink);
            document.write(&mut out)?;
            out.finish()?;
            return Ok(());
        }
        ArchiveEncoding::Text => false,
        ArchiveEncoding::Binary => true,
    };
    // The messages go among the document's bytes, where NCBI's buffer has handed them to its
    // output (`flush`).
    let posts = PostOut {
        diagnostics,
        handle,
    };
    let unposted = write_asn(binary, &mut document, sink, Some(posts))?;
    debug_assert!(unposted.is_empty());
    Ok(())
}

/// The archive of one query batch (`write_search_archive`) as an `ArchiveDocument`.
struct SearchArchive<'h, 'b, A> {
    request: &'h super::ArchiveRequest,
    batch: &'h BatchArchive<'b>,
    search_space_option: Option<i64>,
    database_list: Option<&'h Param>,
    messages: &'h [Message],
    aligns: &'h mut A,
}

impl<A: AlignSource> ArchiveDocument for SearchArchive<'_, '_, A> {
    fn write<E: Encoder>(&mut self, out: &mut ArchiveOut<'_, E>) -> Result<()> {
        write_search_archive(
            out,
            self.request,
            self.batch,
            self.search_space_option,
            self.database_list,
            self.messages,
            &mut *self.aligns,
        )
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
    // The messages are sanitized (`diag_log::sanitize`: 0x20..0x7E only), so no stream writes
    // `#` here: one writing.
    let mut buf = Vec::new();
    let finish = match handle.encoding() {
        ArchiveEncoding::Nothing => return Ok(()),
        ArchiveEncoding::Xml => {
            let mut out = ArchiveOut::once(XmlEncoder::archive(&mut buf), sink);
            write_error_document(&mut out, &messages);
            out.finish()?
        }
        ArchiveEncoding::Text => {
            let mut out = ArchiveOut::once(AsnTextEncoder::new(&mut buf), sink);
            write_error_document(&mut out, &messages);
            out.finish()?
        }
        ArchiveEncoding::Binary => {
            let mut out = ArchiveOut::once(AsnBinaryEncoder::new(&mut buf), sink);
            write_error_document(&mut out, &messages);
            out.finish()?
        }
    };
    debug_assert!(matches!(finish, Finish::Done(ref bad) if bad.is_empty()));
    Ok(())
}
