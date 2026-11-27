//! The annotations of a Seq-entry as NCBI's object manager takes them when a program adds the entry
//! to its scope (line RG D3-1, `om/SPEC.md`; a Bioseq or a Bioseq-set, `entry`): psiblast and
//! tblastn add the query Seq-entry of a PSSM (`index_entry`), makeblastdb each object of its ASN.1
//! input, each in a new scope (`index_object`, line RG D3c-6a). The Seq-tables' column checks when
//! each annotation is added (phase C, `CSeqTableInfo`), the attach of the entry's Bioseq-set and
//! Bioseq ids (`entry::Tse::attach`), then the index of every annotation (phase I,
//! `CSeq_annot_Info::x_InitAnnotKeys`), with NCBI's messages, its exceptions, and the failed
//! diagnostic stream after its "Failed to parse location of" message
//! (`ncbi_exception::main_diag_failbit`), and BitMagic's decoding of a Seq-table's bit vectors
//! (`bitmagic`). What NCBI does not do the same way on every run or every build (a message with a
//! memory address, a crash, BitMagic data that NCBI reads past, decodes with undefined behaviour
//! or from memory it did not write, or answers differently in its scalar and SIMD builds, a signed
//! overflow, the order of memory addresses) is LOSAT's refusal, each naming its condition
//! (`Refusal`).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:535-537 (c++/src/app/blast/tblastn_app.cpp:222-226)
//! NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:494-496
//! NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:230-238 (`CScope_Impl::AddSeq_entry`:
//! the TSE is made, each annotation's `CSeq_annot_Info` with it, then `x_ClearCacheOnNewData`
//! asks for the TSE's annotation Seq-ids, `CTSE_Info::GetAnnotIds`, which indexes them all,
//! `UpdateAnnotIndex`, tse_info.cpp:709, 1035)

mod align;
mod bitmagic;
mod entry;
mod range;
mod table;

use range::{Handle, HandleRangeMap, LiveHandles, Mode, Range};
use view::V;

use crate::ncbi_exception::{post_line_in, NcbiException, Severity};
use crate::report::diag_log::Posts;
use crate::report::serial::asn_types::view;
use crate::report::serial::asn_types::{resolved, table::named, write::failed_stream_text};
use crate::report::serial::pssm::SeqEntry;

/// What NCBI's object manager made of the query's annotations: its messages in order, each one
/// post of its severity (standard error shows their lines; the BLAST archive saves each as one
/// message, `diag_log::Posts::write`), whether its last left the main thread's diagnostic stream
/// failed, and the exception that ends the run (`CATCH_ALL`'s report, already in the failed
/// stream's form when it is).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct AnnotOutcome {
    pub(crate) posts: Posts,
    pub(crate) muted: bool,
    pub(crate) abort: Option<Abort>,
}

/// What ends the run as the query joins NCBI's scope: `CATCH_ALL`'s report and exit status of an
/// exception of the object manager, or LOSAT's refusal of what NCBI does not do the same way on
/// every run (`Refusal`), which stops the run at the same place, after what NCBI writes before
/// (psiblast's prolog of its format; tblastn writes nothing before), in LOSAT's words.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum Abort {
    Native {
        exit: i32,
        message: Vec<u8>,
    },
    /// LOSAT's refusal: its whole text (`in_pssm::PssmReader::refuse`).
    Refused(String),
}

impl Abort {
    /// The end of the run by NCBI's error `error` (its `CATCH_ALL` report), `muted` when the
    /// object manager failed the diagnostic stream: the `LOG_POST` then writes its line feed alone
    /// (`cli::catch_all_text`, which the run applies once the stream is failed; the text kept is
    /// that line feed).
    pub(crate) fn of(error: crate::cli::NativeError, muted: bool) -> Abort {
        Abort::Native {
            exit: error.exit,
            message: if muted { b"\n".to_vec() } else { error.message },
        }
    }

    /// The error that ends the run (`PsiSearch::annot_abort`).
    pub(crate) fn error(&self) -> anyhow::Error {
        match self {
            Abort::Native { exit, message } => crate::cli::NativeError {
                exit: *exit,
                message: message.clone(),
            }
            .into(),
            Abort::Refused(text) => anyhow::anyhow!("{text}"),
        }
    }
}

impl AnnotOutcome {
    /// The outcome of LOSAT's refusal `text` of the query's entry (`Refusal`, `index_entry`):
    /// NCBI's messages before the point it refuses are not written, as its own are not.
    pub(crate) fn refused(text: String) -> AnnotOutcome {
        AnnotOutcome {
            abort: Some(Abort::Refused(text)),
            ..AnnotOutcome::default()
        }
    }
}

/// What LOSAT refuses: an annotation that NCBI's object manager does not index the same way on
/// every run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Refusal {
    /// `Empty region in <address> ...` (1204.1-.4, .7, .8) before the stream failed.
    EmptyRegion,
    /// A null Seq-id handle or location of a table row before the stream failed.
    NullPointer,
    /// A Dense-seg whose division by `numseg` is by zero (SIGFPE).
    DivideByZero,
    /// A bit-bvector or bit-set-bvector that the index reads whose BitMagic data NCBI decodes past
    /// its checks, or that NCBI's builds of BitMagic decode differently or answer differently for
    /// a row of it (`bitmagic::Fail::Refuse`).
    BitVector,
    /// A Packed-seg read past its `present`, `starts`, `ids` or `lens`.
    PackedOutOfBounds,
    /// An int-delta sum that overflows its accumulator.
    DeltaOverflow,
    /// An int-scaled value that overflows `Int8`.
    ScaledOverflow,
    /// A `field-id` whose difference with the column group's base overflows `int`.
    FieldIdOverflow,
    /// Two Seq-ids of one location (a feature's, a table's or a table row's) with different
    /// messages whose order is that of memory addresses (one `m_Packed`).
    IdOrder,
    /// A feature table of more rows than `table::ROWS_LIMIT`.
    TableRows,
    /// A Std-seg `dim` above `align::STD_DIM_LIMIT`.
    StdDim,
    /// An element stated as `NULL` of a list that NCBI's object manager reads through the
    /// reference before the stream failed (`RefElement::deref`), or through a pointer it took
    /// from the list or inside an alignment (`RefElement::deref_raw`).
    NullElement,
    /// A Seq-id handle twice in the query's entry of a search whose scope holds other sequences
    /// first (`-subject`), when a Seq-id of NCBI's message may be written as a subject spelled it
    /// (`entry::Tse::attach`, `seqid::mask_handle::handle_fasta_is_own`).
    DuplicateIdInScope,
}

/// What a refusal's text names (`Refusal::text`), as its caller calls the entry: the entry by its
/// annotations, the entry itself, and the lists with an element stated as `NULL` that the
/// caller's program dereferences (`Refusal::NullElement`: the object manager's, and the caller's
/// own reads of the entry).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Subject {
    /// "a PSSM query whose annot".
    pub(crate) annot: &'static str,
    /// "a PSSM query entry".
    pub(crate) entry: &'static str,
    /// "NCBI's object manager reads when the entry joins its scope (a member of a Bioseq-set,
    /// ...)".
    pub(crate) null_lists: &'static str,
}

impl Refusal {
    /// The refusal's text, with the entry named as `subject` names it, in LOSAT's form (the caller
    /// adds `" is not supported by LOSAT's <PROGRAM>"`: `in_pssm::PssmReader::refuse`, makeblastdb's
    /// `asn_source`).
    pub(crate) fn text(self, subject: Subject) -> String {
        let Subject {
            annot,
            entry,
            null_lists,
        } = subject;
        match self {
            Refusal::EmptyRegion => format!("{annot} has an empty location (NCBI's object manager posts 'Empty region in' with a memory address)"),
            Refusal::NullPointer => format!("{annot} has a Seq-table row without a Seq-id or a location (NCBI's object manager stops with a NULL pointer report of memory addresses)"),
            Refusal::DivideByZero => format!("{annot} has a Dense-seg of no segments with starts or strands (NCBI's object manager divides by zero and crashes)"),
            Refusal::BitVector => format!("{annot} has a Seq-table bit-bvector or bit-set-bvector with malformed BitMagic data that NCBI's object manager reads past, decodes with undefined behaviour (an XOR block, a shift of 32 bits or more, a write past BitMagic's blocks or buffers) or from memory it did not write, or that holds a GAP block that is not one which NCBI's builds of BitMagic (scalar, SSE2, SSE4.2) decode differently or answer differently for a row of the index"),
            Refusal::PackedOutOfBounds => format!("{annot} has a Packed-seg with fewer present, starts, ids or lens than its dim and numseg (NCBI's object manager reads past them)"),
            Refusal::DeltaOverflow => format!("{annot} has a Seq-table int-delta column whose sums leave their integer type (a signed overflow in NCBI's object manager)"),
            Refusal::ScaledOverflow => format!("{annot} has a Seq-table int-scaled column whose values leave Int8 (a signed overflow in NCBI's object manager)"),
            Refusal::FieldIdOverflow => format!("{annot} has a Seq-table column field-id below -2147483638 (a signed overflow in NCBI's object manager)"),
            Refusal::IdOrder => format!("{annot} has a location whose Seq-ids get different messages and have no different gi, accession or general id numbers (NCBI's object manager writes them in the order of memory addresses)"),
            Refusal::TableRows => format!("{annot} has a feature Seq-table of more than 10000000 rows (NCBI's object manager runs out of memory depending on the machine)"),
            Refusal::StdDim => format!("{annot} has a Std-seg dim above 1000000 (NCBI's object manager allocates its rows depending on the machine's memory)"),
            Refusal::NullElement => format!("{entry} with an element stated as NULL in a list that {null_lists}, which it dereferences: a report of memory addresses or a crash"),
            Refusal::DuplicateIdInScope => format!("{entry} with a Seq-id twice in a search of -subject records, when the entry's Bioseqs have a Seq-id that NCBI's duplicate Bioseq id exception may write as a subject record spelled it (a GenBank, EMBL or DDBJ accession alone, a name or a patent with a letter, more than 64 letters in a local or general id or an accession)"),
        }
    }
}

/// What stops the index: an exception of NCBI's (`CException`, a `std::exception`), a null
/// pointer that NCBI's code dereferences through a reference (`CCoreException::eNullPtr`: the
/// refusal, a table row's `Refusal::NullPointer` or a list element's `Refusal::NullElement`,
/// before the stream failed, a report of nothing after it), the failed write of an object's ASN.1 text to the failed stream
/// (`Indexer::empty_region`), or LOSAT's refusal.
#[derive(Debug)]
pub(super) enum Stop {
    Exception(NcbiException),
    Std(&'static str),
    NullPointer(Refusal),
    FailedWrite,
    Refuse(Refusal),
}

impl Stop {
    fn exception(e: NcbiException) -> Stop {
        Stop::Exception(e)
    }

    fn std_exception(what: &'static str) -> Stop {
        Stop::Std(what)
    }
}

/// How NCBI's object manager takes an element of a list of `CRef<>`s of the query's entry: the
/// generic typed access (`asn_types::view::V`, whose `is_null` tells an element stated as `NULL`)
/// with the index's stops for such an element (`Refusal::NullElement`).
pub(super) trait RefElement: Sized {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiobj.cpp:1002-1014
    /// The element of a list of `CRef<>`s that NCBI's code dereferences here through the
    /// reference (`CCoreException::eNullPtr` for one stated as `NULL`: a report of memory
    /// addresses, or of nothing after the stream failed, `Stop::NullPointer`; probe D3c-3b
    /// `probe7`: the parts of a `mix` location, a feature's `ids` and `xref`, an annotation's `id`
    /// and `desc`).
    fn deref(self) -> Result<Self, Stop>;

    /// The element of a list of `CRef<>`s that NCBI's code reads here through a pointer it took
    /// from the list (a crash for one stated as `NULL`: an alignment or a graph of an annotation,
    /// probe D3c-3b `probe7`), or inside an alignment, where LOSAT does not follow which reads
    /// throw and which crash: refused (`Refusal::NullElement`), also after the stream failed.
    fn deref_raw(self) -> Result<Self, Stop>;
}

impl RefElement for V<'_> {
    fn deref(self) -> Result<Self, Stop> {
        if self.is_null() {
            Err(Stop::NullPointer(Refusal::NullElement))
        } else {
            Ok(self)
        }
    }

    fn deref_raw(self) -> Result<Self, Stop> {
        if self.is_null() {
            Err(Stop::Refuse(Refusal::NullElement))
        } else {
            Ok(self)
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:1227-1252 (`x_MapFeatIds`, after
// the feature's locations are mapped: its `id`, each of its `ids`, a gene's names, then each of
// its `xref`s)
/// The elements of the feature `feat` that its index by its ids reads (`x_MapFeatIds`): each of
/// its `ids` and of its `xref`s, through the reference (`RefElement::deref`). The index by ids
/// changes no output of LOSAT's.
fn feature_ids(feat: V<'_>) -> Result<(), Stop> {
    for name in ["ids", "xref"] {
        for element in feat.get(name).into_iter().flat_map(V::items) {
            element.deref()?;
        }
    }
    Ok(())
}

/// The message of an empty overlapping range: the annotation's kind
/// (`seq_annot_info.cpp` `ERR_POST_X(1..8, "Empty region in " ...)`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Event {
    /// 1204.1, a feature (with its ASN.1 text).
    Feature,
    /// 1204.2, a graph.
    Graph,
    /// 1204.3, an alignment.
    Align,
    /// 1204.4, a location set.
    Locs,
    /// 1204.7, a table's location.
    TableLocation,
    /// 1204.8, a feature table's row (no ASN.1 text).
    Row,
}

/// The state of the index: the program's diagnostic prefix, the messages so far and whether the
/// stream has failed.
pub(super) struct Indexer<'p> {
    program: &'p str,
    posts: Posts,
    muted: bool,
    caches: table::Caches,
    /// NCBI's Seq-id handles of the step (`range::LiveHandles`).
    handles: LiveHandles<'p>,
}

/// How the object manager's step ends the run after its messages (`Indexed::end`).
#[derive(Debug)]
pub(crate) enum End {
    /// An exception of NCBI's (`CException`), which `CATCH_ALL` reports.
    Exception(NcbiException),
    /// A `std::exception` and its `what()`.
    Std(&'static str),
    /// A null pointer that NCBI's code dereferences through a reference after the stream failed
    /// (`CCoreException::eNullPtr`): `CATCH_ALL`'s report is its line feed alone (before the
    /// stream failed the report holds memory addresses: `Refusal::NullPointer`,
    /// `Refusal::NullElement`).
    NullPointer,
    /// The failed write of an object's ASN.1 text to the failed stream (`CIOException::eWrite`,
    /// `Indexer::empty_region`), which `CATCH_ALL` ends with nothing more written and the run's
    /// exit status, 0 (`cli::NativeError::silent_end`; `om/` cases f33, k04, k07, k08, k09, k13,
    /// k16, k17; makeblastdb: survey D3c-6 x48, x51).
    FailedWrite,
}

impl End {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-251, 260-263
    /// `CATCH_ALL`'s report and exit status of the end, as the run writes it before the stream
    /// failed (`cli::catch_all_text` turns a report into its line feed after it): of an
    /// exception, of a `std::exception`, the line feed of a null pointer's report after the
    /// stream failed, and no text for the failed write.
    pub(crate) fn native(&self) -> crate::cli::NativeError {
        match self {
            End::Exception(e) => crate::cli::NativeError::exception(e),
            End::Std(what) => crate::cli::NativeError::std_exception(what),
            End::NullPointer => crate::cli::NativeError {
                exit: 255,
                message: b"\n".to_vec(),
            },
            End::FailedWrite => crate::cli::NativeError::silent_end(),
        }
    }
}

/// What NCBI's object manager made of an entry as it joined the scope: its messages in order, each
/// one post of its severity, whether its last left the main thread's diagnostic stream failed, and
/// how it ends the run (`End`), if it does.
#[derive(Debug)]
pub(crate) struct Indexed {
    pub(crate) posts: Posts,
    pub(crate) muted: bool,
    pub(crate) end: Option<End>,
}

/// The object manager's step as the TSE `tse` joins the scope of `program` (the diagnostic
/// prefix): phase C, the attach of its Bioseq-sets' and Bioseqs' ids, phase I, or LOSAT's refusal.
/// `ids_known`: no other sequence is in NCBI's scope before it, so the attach's messages name the
/// entry's own Seq-ids. `muted`: the main thread's diagnostic stream failed before (an earlier
/// object of makeblastdb's input), so its posts are their headers alone from the start.
fn index_tse(
    program: &str,
    tse: &entry::Tse<'_>,
    ids_known: bool,
    muted: bool,
) -> Result<Indexed, Refusal> {
    let mut ix = Indexer {
        program,
        posts: Posts::new(),
        muted,
        caches: table::Caches::default(),
        handles: LiveHandles::default(),
    };
    let end = match ix.run(tse, ids_known) {
        Ok(()) => None,
        Err(Stop::Refuse(refusal)) => return Err(refusal),
        Err(Stop::NullPointer(refusal)) if !ix.muted => return Err(refusal),
        Err(Stop::NullPointer(_)) => Some(End::NullPointer),
        Err(Stop::Exception(e)) => Some(End::Exception(e)),
        Err(Stop::Std(what)) => Some(End::Std(what)),
        Err(Stop::FailedWrite) => Some(End::FailedWrite),
    };
    Ok(Indexed {
        posts: ix.posts,
        muted: ix.muted,
        end,
    })
}

/// The object manager's step as the query's Seq-entry `entry` of a PSSM read by `program`
/// (`psiblast`, `tblastn`: the diagnostic prefix) joins the scope (`index_tse`), or LOSAT's
/// refusal. `ids_known`: no other sequence is in NCBI's scope before it (a search of `-db`), so the
/// attach's messages name the entry's own Seq-ids. `refusal` is LOSAT's refusal of an exception
/// that a check of the data throws where the run's settings turn the check off
/// (`ncbi_environment::ReadVerification`): the run stops there in LOSAT's words, as at a
/// `Refusal`.
pub(crate) fn index_entry(
    program: &str,
    entry: &SeqEntry,
    ids_known: bool,
    refusal: &dyn Fn(&NcbiException) -> Option<String>,
) -> Result<AnnotOutcome, Refusal> {
    let Indexed { posts, muted, end } =
        index_tse(program, &entry::Tse::of(entry), ids_known, false)?;
    let abort = match end {
        None => None,
        Some(End::Exception(e)) => match refusal(&e) {
            Some(text) => return Ok(AnnotOutcome::refused(text)),
            None => Some(Abort::of(crate::cli::NativeError::exception(&e), muted)),
        },
        // The failed write ends the run with no text, whatever the stream's state.
        Some(End::FailedWrite) => {
            let silent = End::FailedWrite.native();
            Some(Abort::Native {
                exit: silent.exit,
                message: silent.message,
            })
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251, 260-263
        // After the stream failed, the `LOG_POST` writes its line feed alone (`cli::catch_all_text`,
        // which the run applies once the stream is failed; the text kept here is that line feed).
        Some(end) => Some(Abort::of(end.native(), muted)),
    };
    Ok(AnnotOutcome {
        posts,
        muted,
        abort,
    })
}

/// The object manager's step as an object of makeblastdb's ASN.1 input joins a new scope of its
/// own (`index_tse`; `CSeqEntrySource`): `object` is the Seq-entry, or the Bioseq-set when
/// `is_set`, as the reader read it, and `muted` whether an earlier object left the main thread's
/// diagnostic stream failed. No other sequence is in the scope.
pub(crate) fn index_object(object: V<'_>, is_set: bool, muted: bool) -> Result<Indexed, Refusal> {
    index_tse(
        "makeblastdb",
        &entry::Tse::of_object(object, is_set),
        true,
        muted,
    )
}

/// `index_entry` of a query Bioseq with the Seq-ids `query_ids` and the annotations `annots`.
#[cfg(test)]
pub(crate) fn index_annots(
    program: &str,
    annots: &[crate::report::serial::asn_types::Value],
    query_ids: &[crate::seqid::ParsedSeqId],
) -> Result<AnnotOutcome, Refusal> {
    use crate::report::serial::pssm::{PssmQuery, SeqInst};
    let query = PssmQuery {
        ids: query_ids.to_vec(),
        descr: None,
        inst: SeqInst::new(("raw", 2), ("aa", 3), None, None),
        annot: Some(annots.to_vec()),
    };
    index_entry(program, &SeqEntry::Seq(Box::new(query)), true, &|_| None)
}

impl<'p> Indexer<'p> {
    /// A post of the object manager (`ERR_POST_X`, `ERR_POST`; `id` empty without an error code),
    /// in the stream's state.
    fn post(&mut self, severity: Severity, id: &str, text: &[u8]) {
        let id = (!id.is_empty()).then_some(id);
        let code = match severity {
            Severity::Warning => crate::report::diag_log::WARNING,
            Severity::Error => crate::report::diag_log::ERROR,
        };
        self.posts.post(
            code,
            post_line_in(severity, id, Some(self.program), text, self.muted),
        );
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:403-427,474-500,646-663 (phase C:
    // `x_InitAnnotList` of each annotation in order makes a Seq-table's `CSeqTableInfo`; the other
    // kinds only list their objects, a feature's by its type, through the reference) and :639-668
    // (phase I, `x_InitAnnotKeys`)
    // NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:98-111 (a feature's type, read
    // when phase C lists it)
    fn run(&mut self, tse: &entry::Tse<'p>, ids_known: bool) -> Result<(), Stop> {
        let annots = &tse.annots;
        let data: Vec<Option<(&'static str, V<'_>)>> = annots
            .iter()
            .map(|annot| annot.get("data").map(V::choice))
            .collect();
        // A table's information is boxed: one slot of each annotation is a pointer, not the
        // size of a `TableInfo` (hundreds of bytes per annotation of a feature-rich entry).
        let mut tables: Vec<Option<Box<table::TableInfo<'_>>>> =
            std::iter::repeat_with(|| None).take(annots.len()).collect();
        let mut ids: Vec<Vec<Handle>> = vec![Vec::new(); tse.objects()];
        for &made in &tse.construction {
            match made {
                entry::Made::Bioseq(object) => {
                    ids[object] = tse
                        .bioseq_ids(object)
                        .iter()
                        .map(|id| Handle::of_run(id, &mut self.handles))
                        .collect();
                }
                entry::Made::Annot(index) => match data[index] {
                    Some(("seq-table", table)) => {
                        tables[index] = Some(Box::new(self.table_info(table)?));
                    }
                    Some(("ftable", feats)) => {
                        for feat in feats.items() {
                            feat.deref()?;
                        }
                    }
                    _ => {}
                },
                entry::Made::Null => return Err(Stop::NullPointer(Refusal::NullElement)),
            }
        }
        tse.attach(&ids, &self.handles, ids_known)?;
        for ((annot, kind), table) in annots.iter().zip(&data).zip(&tables) {
            match kind {
                Some(("ftable", feats)) => {
                    for feat in feats.items() {
                        self.index_feature(feat)?;
                    }
                }
                Some(("align", aligns)) => {
                    for align in aligns.items() {
                        self.index_align(align.deref_raw()?)?;
                    }
                }
                Some(("graph", graphs)) => {
                    for graph in graphs.items() {
                        let graph = graph.deref_raw()?;
                        let loc = graph.field("loc");
                        let map = HandleRangeMap::of_location(loc, Mode::None, &mut self.handles)?;
                        self.refuse_empty(std::slice::from_ref(&map), Event::Graph, graph)?;
                        self.handles.end_object();
                    }
                }
                Some(("locs", locs)) => self.index_locs(*annot, *locs)?,
                Some(("seq-table", _)) => self.index_table(table.as_deref().expect("phase C"))?,
                // `ids` (and an annotation without data) are not indexed.
                _ => {}
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:426-447 (`x_ProcessFeat`)
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:741-786 (`x_InitFeatKeys`)
    fn index_feature(&mut self, feat: V<'_>) -> Result<(), Stop> {
        let except = feat.get("except-text").map(V::str).unwrap_or_default();
        let mode = if contains(except, b"trans-splicing") {
            Mode::TransSplicing
        } else if contains(except, b"circular RNA") {
            Mode::CircularRna
        } else {
            Mode::None
        };
        // Both maps are made (`GetMaps`, the dereferences of their parts) before either is mapped.
        let mut maps: [Option<HandleRangeMap>; 2] = [None, None];
        for (map, loc) in maps
            .iter_mut()
            .zip([feat.get("location"), feat.get("product")])
        {
            if let Some(loc) = loc {
                *map = Some(HandleRangeMap::of_location(loc, mode, &mut self.handles)?);
            }
        }
        for map in maps.iter().flatten() {
            self.map_location(map, Event::Feature, feat)?;
        }
        feature_ids(feat)?;
        self.handles.end_object();
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:845-915 (`x_AddAlignKeys`,
    // `x_AddAlignIndex`)
    // The mapped range is never empty here (no "Failed to parse location"): an alignment posts
    // only its sizes' warnings, or its empty regions (refused), in any order of its Seq-ids. The
    // Seq-ids of one row, and for a discontinuous alignment of all its rows, are one handle as
    // NCBI's trees take them (`Handle`).
    fn index_align(&mut self, align: V<'_>) -> Result<(), Stop> {
        let segs = align.get("segs").map(V::choice);
        if let Some(("disc", subs)) = segs {
            let mut all = Vec::new();
            for sub in subs.items() {
                let mut rows = Vec::new();
                align::process_align(self, &mut rows, sub.deref_raw()?)?;
                all.push(rows);
            }
            let mut total: Vec<(&Handle, Range)> = Vec::new();
            for map in all.iter().flatten() {
                for (handle, hr) in &map.map {
                    let range = hr.overlapping();
                    match total.iter_mut().find(|(h, _)| *h == handle) {
                        Some((_, sum)) => sum.combine(range),
                        None => total.push((handle, range)),
                    }
                }
            }
            if total.iter().any(|(_, range)| range.is_empty()) {
                return Err(self.empty_region(Event::Align, align));
            }
            self.handles.mapped_all(all.iter().flatten());
        } else {
            let mut rows = Vec::new();
            align::process_align(self, &mut rows, align)?;
            for map in &rows {
                self.refuse_empty(std::slice::from_ref(map), Event::Align, align)?;
            }
        }
        self.handles.end_object();
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:940-975 (`x_InitLocsKeys`)
    // NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:396-416 (`GetMaps` of `locs`: the
    // first `region` of the annotation's `desc`, not the location itself)
    fn index_locs(&mut self, annot: V<'_>, locs: V<'_>) -> Result<(), Stop> {
        if locs.len() != 1 {
            return Ok(());
        }
        let region = annot.get("desc").and_then(|desc| {
            desc.items()
                .map(V::choice)
                .find(|(variant, _)| *variant == "region")
                .map(|(_, region)| region)
        });
        if let Some(region) = region {
            let map = HandleRangeMap::of_location(region, Mode::None, &mut self.handles)?;
            // The message writes the location set's one location (`info.GetLocs()`), not the
            // region.
            let loc = locs.item(0).expect("one location");
            self.refuse_empty(std::slice::from_ref(&map), Event::Locs, loc)?;
            self.handles.end_object();
        }
        Ok(())
    }

    /// An empty region of an annotation whose index maps only non-empty ranges (a graph, an
    /// alignment's row, a location set; `event` names its message, `object` what it writes):
    /// the same in any order of the Seq-ids (`empty_region`).
    fn refuse_empty(
        &mut self,
        maps: &[HandleRangeMap],
        event: Event,
        object: V<'_>,
    ) -> Result<(), Stop> {
        debug_assert!(matches!(event, Event::Graph | Event::Align | Event::Locs));
        let empty =
            |map: &HandleRangeMap| map.map.iter().any(|(_, hr)| hr.overlapping().is_empty());
        if maps.iter().any(empty) {
            return Err(self.empty_region(event, object));
        }
        self.handles.mapped_all(maps);
        Ok(())
    }

    /// `ERR_POST_X` (subcode `n`) of "Empty region in ", `GetDescription()`, a space and the
    /// ASN.1 text (`MSerial_AsnText`) of `object`, for an empty region of `object` (`event`, not a
    /// table row: `map_location`): before the stream failed, a message with the TSE's memory
    /// address (`CTSE_Info::GetDescription`, tse_info.cpp:2049-2056), refused. After it, the post's
    /// header and what the failed stream gets of the object's ASN.1 text
    /// (`asn_types::write::failed_stream_text`: the text without its last end of line, or its start
    /// where NCBI's buffer of 4096 bytes first flushed), whose write throws (`Stop::FailedWrite`;
    /// `om/SPEC.md` §3 rule 3).
    fn empty_region(&mut self, event: Event, object: V<'_>) -> Stop {
        if !self.muted {
            return Stop::Refuse(Refusal::EmptyRegion);
        }
        // The one location of a `locs` stated as `NULL`, which the message writes.
        if object.is_null() {
            return Stop::Refuse(Refusal::NullElement);
        }
        let (id, root) = match event {
            Event::Feature => ("1204.1", named::SEQ_FEAT),
            Event::Graph => ("1204.2", named::SEQ_GRAPH),
            Event::Align => ("1204.3", named::SEQ_ALIGN),
            Event::Locs => ("1204.4", named::SEQ_LOC),
            Event::TableLocation => ("1204.7", named::SEQ_LOC),
            Event::Row => unreachable!("a table row's empty region writes no object"),
        };
        debug_assert!(std::ptr::eq(resolved(object.ty), resolved(root)));
        let mut line = post_line_in(Severity::Error, Some(id), Some(self.program), b"", true);
        line.pop();
        line.extend_from_slice(&failed_stream_text(root, object.value));
        line.push(b'\n');
        self.posts.post(crate::report::diag_log::ERROR, line);
        Stop::FailedWrite
    }

    /// The null pointer of a table row (`Stop::NullPointer`).
    fn null_pointer(&self) -> Stop {
        Stop::NullPointer(Refusal::NullPointer)
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:741-786 (cited at
    // `index_feature`), 1048-1080 and 1119-1150 (a table's location and a table row's location:
    // the same loop with `ERR_POST_X` of subcode 7, and of subcode 8 that posts "Empty region in ",
    // `GetDescription()`, " Seq-table row " and the row)
    // Each Seq-id of the map in NCBI's order (`HandleRangeMap::ncbi_order`): its overlapping range,
    // or for a circular one its two ranges, mapped (`x_map`); an empty one posts "Empty region"
    // (`empty_region`: refused before the stream failed, a memory address; after it the object's
    // ASN.1 text, which ends the run), except for a table row after the stream failed, whose post
    // is its header alone, `Error: (1204.8)`. Where NCBI's order is that of memory addresses
    // (Seq-ids of one `m_Packed`: by `m_Info`, the lower address first,
    // seq_id_handle.hpp:207, cited at `HandleRangeMap::ncbi_order`), the order changes the messages
    // only when an empty Seq-id's post and another's "Failed to parse location" (`Error: (1204.6)`)
    // are in one group: refused (`om/` cases `s_o09_tie_gi5_acc5_u` and
    // `s_o10_tie_acc5v1_acc5v2_u`: psiblast and tblastn wrote a row's in different orders; in
    // another kind of annotation the 1204.6 message comes before the object's text, or the run
    // stops on the empty region before it). Any other messages of a group are the same in every
    // order (the 1204.6 messages differ only by how many came before, the object's text is that of
    // the object, not of a Seq-id).
    fn map_location(
        &mut self,
        map: &HandleRangeMap,
        event: Event,
        object: V<'_>,
    ) -> Result<(), Stop> {
        let empty = |hr: &range::HandleRange| hr.overlapping().is_empty();
        let fails = |hr: &range::HandleRange| {
            !empty(hr)
                && hr.has_gaps()
                && hr.circular
                && (hr.circular_start().is_empty() || hr.circular_end().is_empty())
        };
        for group in map.ncbi_order() {
            if group.iter().any(|(_, hr)| empty(hr)) && group.iter().any(|(_, hr)| fails(hr)) {
                return Err(Stop::Refuse(Refusal::IdOrder));
            }
            for (handle, hr) in group {
                let range = hr.overlapping();
                if range.is_empty() {
                    if event != Event::Row {
                        return Err(self.empty_region(event, object));
                    }
                    if !self.muted {
                        return Err(Stop::Refuse(Refusal::EmptyRegion));
                    }
                    // After the stream failed: the header alone.
                    self.post(Severity::Error, "1204.8", b"");
                    continue;
                }
                if hr.has_gaps() && hr.circular {
                    self.x_map(hr.circular_start(), Some(handle.clone()))?;
                    self.x_map(hr.circular_end(), Some(handle.clone()))?;
                } else {
                    self.x_map(range, Some(handle.clone()))?;
                }
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:672-700
    // The message stops at `s.rdbuf()`, which fails the stream (`ncbi_exception::
    // main_diag_failbit`). `mapper.Map` of a null handle (a table row without a Seq-id) asks
    // `CSeq_id_Handle::HaveMatchingHandles()` of it: `CCoreException::eNullPtr` (`Stop::NullPointer`).
    /// `x_Map` of `range` and `handle`: `None` is a null handle.
    fn x_map(&mut self, range: Range, handle: Option<Handle>) -> Result<(), Stop> {
        if range.is_empty() {
            self.post(Severity::Error, "1204.6", b"Failed to parse location of ");
            self.muted = true;
            return Ok(());
        }
        let Some(handle) = handle else {
            return Err(self.null_pointer());
        };
        self.handles.mapped(&handle);
        Ok(())
    }
}

/// `string::find(needle) != NPOS`.
fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    haystack
        .windows(needle.len())
        .any(|window| window == needle)
}

#[cfg(test)]
mod tests;
