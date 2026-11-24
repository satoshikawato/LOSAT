//! The reader of makeblastdb's mask files (`makeblastdb.cpp` `s_ReadObject` and
//! `x_ProcessMaskData`): the format guess, then `Blast-db-mask-info` objects one after another
//! in ASN.1 text or binary.
//!
//! NCBI reads an object with its serial library (`CObjectIStreamAsn`, `CObjectIStreamAsnBinary`,
//! driven by the type information of the generated classes). LOSAT reads it with the reader of
//! the type table (`report::serial::asn_types`: the object stack of `asn_types::read`, NCBI's
//! frames, members, pointer forms and errors, over the primitives of either encoding), whose
//! typed calls (`asn_types::typed`) build makeblastdb's model of the object as it is read
//! (`Build`; the Seq-id by `asn_types::seq_id`), object by object from the file's stream.
//!
//! What NCBI's reader throws for a file it cannot read is what LOSAT reports (`ReadError`, then
//! `MaskFileError::Ncbi`): the exception of the primitive at its place in NCBI's source, with
//! what each object frame it leaves adds, as makeblastdb's `CATCH_ALL` prints it. What it posts
//! while it reads (a byte of a VisibleString read as `#`, a UTF8String tag) is handed out in order
//! (`MaskFileReader::take_posts`).
//!
//! The end of the file is NCBI's `CEofException` (`ReadError::Eof`). Every object frame turns it
//! into a `CSerialException` (`HandleEOF`) unless fewer than two frames remain after it, in which
//! case the `CEofException` goes on to `x_ProcessMaskData`, which treats it as the end of the file
//! and drops the object being read. That happens at the type name of a text object, at its opening
//! brace, and at the tag, the length and the closing bytes of the outermost binary SEQUENCE.

use std::fmt;
use std::io::Read;
use std::path::{Path, PathBuf};

use super::guess::guess_format;
use super::{MaskInfo, MaskInterval, MaskList, MaskLoc, MaskPart, NullElement};
use crate::report::serial::asn_ber_reader::BerReader;
use crate::report::serial::asn_input::{Input, Names, TextInput, R};
use crate::report::serial::asn_text_reader::StreamSource;
use crate::report::serial::asn_types::read::{ObjectStream, Post, ReadError, TypedReader};
#[cfg(test)]
use crate::report::serial::asn_types::seq_id::read_seq_id;
use crate::report::serial::asn_types::seq_id::{mandatory, read_seq_id_into};
use crate::report::serial::asn_types::table::named;
use crate::seqid::model::ParsedSeqId;
use crate::writedb::WriteDbError;

/// The encoding of a mask file, as `CFormatGuess` finds it. Anything else is an error when the
/// file is opened ("Unknown encoding").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MaskFileFormat {
    /// `eBinaryASN`: `MSerial_AsnBinary`.
    BinaryAsn,
    /// `eTextASN`: `MSerial_AsnText`.
    TextAsn,
}

/// Why a mask file could not be read.
#[derive(Debug)]
pub(crate) enum MaskFileError {
    /// What makeblastdb ends with, as NCBI's code writes it: the "Unknown encoding ..."
    /// `CInvalidDataException` of the format guess, and what NCBI's serial reader throws for a
    /// file it cannot read: its `CSerialException` or `CUtilException` (`Exception`, "Error: "
    /// and the report, status 255), `std::bad_alloc` of a string's length (`OutOfMemory`,
    /// status 4) or `std::length_error` (`Other` with its `what()`, status 255).
    Ncbi(WriteDbError),
    /// NCBI reads it; LOSAT does not.
    Unsupported(String),
    /// The file ends where `next_list` expected the continuation that `more TRUE` announces.
    /// `x_ProcessMaskData` has no handler for the `CEofException` that this raises outside
    /// of its first `try` block: `CATCH_ALL` ignores a `CIOException` that is not a flush
    /// failure, so makeblastdb leaves with status 0 and no database, after the
    /// "No sequences matched any of the masks provided." message of its build object.
    EofInContinuation { path: PathBuf },
}

impl fmt::Display for MaskFileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MaskFileError::Ncbi(e) => write!(f, "{e}"),
            MaskFileError::Unsupported(what) => {
                write!(f, "mask data that LOSAT does not read: {what}")
            }
            MaskFileError::EofInContinuation { path } => write!(
                f,
                "mask data in '{}': the file ends where the continuation of 'more TRUE' was expected",
                path.display()
            ),
        }
    }
}

impl std::error::Error for MaskFileError {}

/// How far `CSeq_loc::CheckId` walks one Seq-loc (`x_CheckId`, cited at `read_seq_loc_parts`):
/// each Seq-loc starts a walk of its own (`my_id`, NULL), takes the Seq-ids of its parts in
/// order, and hands its own to the Seq-loc around it (`x_UpdateId` of `id` and `my_id`) once its
/// walk is done. makeblastdb's `GetId` of a mask is this walk of the mask's Seq-loc.
#[derive(Debug)]
enum IdWalk {
    /// Every Seq-id met so far is `Equals` to this one (`None`: none yet).
    Same(Option<ParsedSeqId>),
    /// `x_CheckId` returned false: two Seq-ids differ, or a `feat` came. Nothing after it is
    /// walked, and `GetId` returns NULL.
    Stopped,
    /// The walk reached an element that is an empty `CRef<>`: `CCoreException` `eNullPtr`.
    Null(NullElement),
}

impl IdWalk {
    /// `x_UpdateId` of `id` and the address of the Seq-id in `seq_id`, for a part: the first
    /// one is taken from the slot, a later one is compared where it is.
    fn update(&mut self, seq_id: &mut Option<ParsedSeqId>) {
        match self {
            IdWalk::Same(None) => *self = IdWalk::Same(seq_id.take()),
            // `Equals` compares the values: the typed Seq-ids are equal exactly when it does.
            IdWalk::Same(Some(first)) if seq_id.as_ref() != Some(first) => *self = IdWalk::Stopped,
            _ => {}
        }
    }

    /// `CheckId` of `id` and `may_throw`, called on an element of a `mix` or `equiv`:
    /// the element's own walk `inner`, if this walk gets that far.
    fn nested(&mut self, inner: IdWalk) {
        if !matches!(self, IdWalk::Same(_)) {
            return;
        }
        match inner {
            IdWalk::Same(Some(seq_id)) => self.update(&mut Some(seq_id)),
            // `x_UpdateId` of `id` and `NULL` keeps `id`.
            IdWalk::Same(None) => {}
            stopped_or_null => *self = stopped_or_null,
        }
    }

    /// An element that is an empty `CRef<>`, dereferenced if this walk gets that far.
    fn null(&mut self, null: NullElement) {
        if matches!(self, IdWalk::Same(_)) {
            *self = IdWalk::Null(null);
        }
    }
}

/// The lists of `CRef`s that a mask file holds, as `NullElement::list` names them.
const LIST_MASKS: &str = "the masks of a Blast-mask-list";
const LIST_MIX: &str = "a Seq-loc-mix (mix)";
const LIST_EQUIV: &str = "a Seq-loc-equiv (equiv)";
const LIST_PACKED_INT: &str = "a Packed-seqint (packed-int)";

/// The read of a mask file's object into makeblastdb's model, on the object stack of the table's
/// reader (`asn_types::typed`): the members, variants and elements that makeblastdb uses are
/// built here, the frames, members, pointer forms and errors are the reader's. The Seq-id of a
/// part is read into a slot of its Seq-loc (`read_seq_id_into`), where the walk takes or compares
/// it (`IdWalk::update`); a part's Seq-id is a mandatory member, so the slot holds it when the
/// part's class was read. The closures of the parts that a mask file repeats (a Seq-loc, an
/// interval or point, an element of `packed-int`) are inlined into the typed calls, as the
/// functions of the reader before them were.
struct Build;

impl Build {
    /// A `Seq-interval` (of the type `ty`), its Seq-id in `id`.
    fn seq_interval<I: Input>(
        r: &mut TypedReader<I>,
        ty: u16,
        id: &mut Option<ParsedSeqId>,
    ) -> R<MaskInterval> {
        let (mut from, mut to, mut strand) = (0u32, 0u32, None);
        *id = None;
        r.class(
            ty,
            #[inline(always)]
            |r, index, ty| {
                match index {
                    0 => from = r.integer(ty)? as u32,
                    1 => to = r.integer(ty)? as u32,
                    // The values of Na-strand are 0..=4 and 255.
                    2 => strand = Some(r.enumerated(ty)? as u8),
                    3 => read_seq_id_into(r, id)?,
                    // fuzz-from, fuzz-to: read, not used.
                    _ => drop(r.value(ty)?),
                }
                Ok(())
            },
        )?;
        debug_assert!(
            id.is_some(),
            "ReadClassSequential has read every mandatory member"
        );
        Ok(MaskInterval { from, to, strand })
    }

    /// A `Seq-point` (also the members `a` and `b` of a `Seq-bond`), its Seq-id in `id`.
    fn seq_point<I: Input>(
        r: &mut TypedReader<I>,
        ty: u16,
        id: &mut Option<ParsedSeqId>,
    ) -> R<MaskInterval> {
        let (mut point, mut strand) = (0u32, None);
        *id = None;
        r.class(
            ty,
            #[inline(always)]
            |r, index, ty| {
                match index {
                    0 => point = r.integer(ty)? as u32,
                    1 => strand = Some(r.enumerated(ty)? as u8),
                    2 => read_seq_id_into(r, id)?,
                    // fuzz: read, not used.
                    _ => drop(r.value(ty)?),
                }
                Ok(())
            },
        )?;
        debug_assert!(
            id.is_some(),
            "ReadClassSequential has read every mandatory member"
        );
        Ok(MaskInterval {
            from: point,
            to: point,
            strand,
        })
    }

    /// A `Packed-seqpnt`: its points, its Seq-id in `seq_id`.
    fn packed_seqpnt<I: Input>(
        r: &mut TypedReader<I>,
        ty: u16,
        parts: &mut Vec<MaskPart>,
        seq_id: &mut Option<ParsedSeqId>,
    ) -> R<()> {
        let (mut strand, mut points) = (None, Vec::new());
        *seq_id = None;
        r.class(ty, |r, member, ty| {
            match member {
                0 => strand = Some(r.enumerated(ty)? as u8),
                1 => read_seq_id_into(r, seq_id)?,
                3 => r.list(ty, |r, element| {
                    points.push(r.integer(element.expect("a list of numbers"))? as u32);
                    Ok(())
                })?,
                // fuzz: read, not used.
                _ => drop(r.value(ty)?),
            }
            Ok(())
        })?;
        parts.extend(points.into_iter().map(|p| {
            MaskPart::Interval(MaskInterval {
                from: p,
                to: p,
                strand,
            })
        }));
        debug_assert!(
            seq_id.is_some(),
            "ReadClassSequential has read every mandatory member"
        );
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:3528-3620
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:403-422
    // NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_loc.hpp:926-945
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:1439-1555,1559-1600
    /// One `Seq-loc` (any variant, nested `mix` and `equiv` included): appends the ranges
    /// `CSeq_loc_CI` visits to `parts`, in its order (a `whole` Seq-loc is the interval 0 to
    /// 4294967294, whose `CRange` is `TRange::GetWhole()`; a point the interval of one position),
    /// and walks it in `walk` (a new [`IdWalk`]) as `CSeq_loc::CheckId` does. An element of
    /// `packed-int`, `mix` or `equiv` stated as `NULL` (an empty `CRef<>`) is where the walk
    /// dereferences it, if it gets that far (through the element's reference, for its range or
    /// for its own `CheckId`).
    fn seq_loc_parts<I: Input>(
        r: &mut TypedReader<I>,
        parts: &mut Vec<MaskPart>,
        walk: &mut IdWalk,
    ) -> R<()> {
        let mut seq_id = None;
        r.choice(
            named::SEQ_LOC,
            #[inline(always)]
            |r, index, ty| {
                match index {
                    SEQ_LOC_NULL => {
                        r.null()?;
                        parts.push(MaskPart::Null);
                    }
                    SEQ_LOC_EMPTY => {
                        read_seq_id_into(r, &mut seq_id)?;
                        walk.update(&mut seq_id);
                        parts.push(MaskPart::Empty);
                    }
                    SEQ_LOC_WHOLE => {
                        read_seq_id_into(r, &mut seq_id)?;
                        walk.update(&mut seq_id);
                        parts.push(MaskPart::Interval(MaskInterval {
                            from: 0,
                            to: u32::MAX - 1,
                            strand: None,
                        }));
                    }
                    SEQ_LOC_INT => {
                        let interval = Self::seq_interval(r, ty, &mut seq_id)?;
                        walk.update(&mut seq_id);
                        parts.push(MaskPart::Interval(interval));
                    }
                    SEQ_LOC_PACKED_INT => r.list(
                        ty,
                        #[inline(always)]
                        |r, element| {
                            let Some(element) = element else {
                                walk.null(null_element(r, LIST_PACKED_INT));
                                return Ok(());
                            };
                            let interval = Self::seq_interval(r, element, &mut seq_id)?;
                            walk.update(&mut seq_id);
                            parts.push(MaskPart::Interval(interval));
                            Ok(())
                        },
                    )?,
                    SEQ_LOC_PNT => {
                        let interval = Self::seq_point(r, ty, &mut seq_id)?;
                        walk.update(&mut seq_id);
                        parts.push(MaskPart::Interval(interval));
                    }
                    SEQ_LOC_PACKED_PNT => {
                        Self::packed_seqpnt(r, ty, parts, &mut seq_id)?;
                        walk.update(&mut seq_id);
                    }
                    SEQ_LOC_MIX | SEQ_LOC_EQUIV => {
                        let list = if index == SEQ_LOC_MIX {
                            LIST_MIX
                        } else {
                            LIST_EQUIV
                        };
                        r.list(ty, |r, element| {
                            if element.is_none() {
                                walk.null(null_element(r, list));
                                return Ok(());
                            }
                            let mut inner = IdWalk::Same(None);
                            Self::seq_loc_parts(r, parts, &mut inner)?;
                            walk.nested(inner);
                            Ok(())
                        })?;
                    }
                    SEQ_LOC_BOND => {
                        // `b` is an OPTIONAL `CRef<>` (`TypedReader::class`): a `NULL` leaves it
                        // unset, which `x_CheckId` and `CSeq_loc_CI` skip (`CanGetB`, `IsSetB`).
                        r.class(ty, |r, _, ty| {
                            let interval = Self::seq_point(r, ty, &mut seq_id)?;
                            walk.update(&mut seq_id);
                            parts.push(MaskPart::Interval(interval));
                            Ok(())
                        })?;
                    }
                    _ => {
                        // `feat`: read; `x_CheckId` fails, so the Seq-loc has no Seq-id.
                        drop(r.value(ty)?);
                        *walk = IdWalk::Stopped;
                    }
                }
                Ok(())
            },
        )
    }

    /// One `Seq-loc` of `masks`: the ranges `CSeq_loc_CI` visits and `CSeq_loc::GetId()`, or the
    /// element stated as `NULL` that `GetId` dereferences.
    fn seq_loc<I: Input>(r: &mut TypedReader<I>) -> R<Result<MaskLoc, NullElement>> {
        let mut parts = Vec::new();
        let mut walk = IdWalk::Same(None);
        Self::seq_loc_parts(r, &mut parts, &mut walk)?;
        let id = match walk {
            IdWalk::Same(id) => id,
            IdWalk::Stopped => None,
            IdWalk::Null(null) => return Ok(Err(null)),
        };
        Ok(Ok(MaskLoc { id, parts }))
    }

    /// A `Blast-mask-list` (of the type `ty`). Its masks are kept up to the first one whose
    /// `GetId` dereferences an element stated as `NULL` (the mask itself, or one inside it), which
    /// is `MaskList::null_element`; the rest is read and dropped.
    fn mask_list<I: Input>(r: &mut TypedReader<I>, ty: u16) -> R<MaskList> {
        let mut masks = Vec::new();
        let mut null_element_at = None;
        let mut more = None;
        r.class(ty, |r, index, ty| {
            if index == 0 {
                return r.list(ty, |r, element| {
                    if element.is_none() {
                        if null_element_at.is_none() {
                            null_element_at = Some(null_element(r, LIST_MASKS));
                        }
                        return Ok(());
                    }
                    let loc = Self::seq_loc(r)?;
                    if null_element_at.is_none() {
                        match loc {
                            Ok(loc) => masks.push(loc),
                            Err(null) => null_element_at = Some(null),
                        }
                    }
                    Ok(())
                });
            }
            more = Some(r.boolean()?);
            Ok(())
        })?;
        Ok(MaskList {
            masks,
            more: mandatory(more),
            null_element: null_element_at,
        })
    }

    /// A `Blast-db-mask-info`.
    fn info<I: Input>(r: &mut TypedReader<I>) -> R<MaskInfo> {
        let mut algo_id = 0;
        let mut algo_program = 0;
        let mut algo_options = Vec::new();
        let mut masks = None;
        r.class(named::BLAST_DB_MASK_INFO, |r, index, ty| {
            match index {
                0 => algo_id = r.integer(ty)? as i32,
                1 => algo_program = r.enumerated(ty)?,
                2 => algo_options = r.text(ty)?,
                _ => masks = Some(Self::mask_list(r, ty)?),
            }
            Ok(())
        })?;
        Ok(MaskInfo {
            algo_id,
            algo_program,
            algo_options,
            masks: mandatory(masks),
        })
    }

    /// A `Blast-mask-list` read as the root (the continuation).
    fn list<I: Input>(r: &mut TypedReader<I>) -> R<MaskList> {
        Self::mask_list(r, named::BLAST_MASK_LIST)
    }
}

/// The element stated as `NULL` that the reader `r` just read, of the list `list`.
fn null_element<I: Input>(r: &TypedReader<I>, list: &'static str) -> NullElement {
    NullElement {
        list,
        position: r.position().to_string(),
    }
}

/// The variants of `Seq-loc` (seqloc.asn), checked against the table (`tests`).
const SEQ_LOC_NULL: usize = 0;
const SEQ_LOC_EMPTY: usize = 1;
const SEQ_LOC_WHOLE: usize = 2;
const SEQ_LOC_INT: usize = 3;
const SEQ_LOC_PACKED_INT: usize = 4;
const SEQ_LOC_PNT: usize = 5;
const SEQ_LOC_PACKED_PNT: usize = 6;
const SEQ_LOC_MIX: usize = 7;
const SEQ_LOC_EQUIV: usize = 8;
const SEQ_LOC_BOND: usize = 9;

/// The member names of the mask types in the order `Build` takes them: the table's module order
/// (`tests`).
#[cfg(test)]
pub(super) const MEMBER_ORDER: [(u16, &[&str]); 7] = [
    (
        named::SEQ_LOC,
        &[
            "null",
            "empty",
            "whole",
            "int",
            "packed-int",
            "pnt",
            "packed-pnt",
            "mix",
            "equiv",
            "bond",
            "feat",
        ],
    ),
    (
        named::SEQ_INTERVAL,
        &["from", "to", "strand", "id", "fuzz-from", "fuzz-to"],
    ),
    (named::SEQ_POINT, &["point", "strand", "id", "fuzz"]),
    (named::PACKED_SEQPNT, &["strand", "id", "fuzz", "points"]),
    (named::SEQ_BOND, &["a", "b"]),
    (named::BLAST_MASK_LIST, &["masks", "more"]),
    (
        named::BLAST_DB_MASK_INFO,
        &["algo-id", "algo-program", "algo-options", "masks"],
    ),
];

/// One `Seq-id` of an ASN.1 binary value, read by the reader of the mask files (`read_seq_id`
/// on the object stack), for the test that compares it with `seqid::ber::decode_seq_id`
/// (DECISIONS 35).
#[cfg(test)]
pub(super) fn seq_id_from_binary(data: &[u8]) -> Result<ParsedSeqId, ReadError> {
    use crate::report::serial::asn_text_reader::SliceSource;
    let mut stream = ObjectStream::binary(SliceSource::new(data), Names::MASK);
    stream.read_typed(named::SEQ_ID, read_seq_id).0
}

/// [`seq_id_from_binary`] for the value `data` in ASN.1 text, for the test that compares the
/// reader with `seqid::asn_text`.
#[cfg(test)]
pub(super) fn seq_id_from_text(data: &[u8]) -> Result<ParsedSeqId, ReadError> {
    use crate::report::serial::asn_text_reader::SliceSource;
    let document = [b"Seq-id ::= ".as_slice(), data].concat();
    let mut stream = ObjectStream::text(SliceSource::new(&document), Names::MASK);
    stream.read_typed(named::SEQ_ID, read_seq_id).0
}

/// Where the reader is in its file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    Open,
    /// The clean end of the file.
    Ended,
    /// An error: every later call is an error again, never a clean end.
    Failed,
    /// The file ended where a continuation was expected.
    FailedEof,
}

fn unknown_encoding(name: &Path) -> MaskFileError {
    MaskFileError::Ncbi(WriteDbError::Options(format!(
        "Unknown encoding for mask data in '{}'",
        name.display()
    )))
}

/// The bytes of a mask file, read in blocks as its objects are read.
type FileSource = StreamSource<Box<dyn Read>>;

/// Where the objects of a mask file come from: one stream for all objects of the file (only the
/// bytes of the object being read are kept), each object read by a new state of NCBI's stream
/// (`ObjectStream`).
enum Source {
    /// ASN.1 text.
    Text(Box<ObjectStream<TextInput<FileSource>>>),
    /// ASN.1 binary.
    Binary(Box<ObjectStream<BerReader<FileSource>>>),
}

/// A mask file read from the disk.
pub(crate) struct MaskFileReader {
    path: PathBuf,
    source: Source,
    /// Set by the first error or end of file (NCBI stops reading the file there).
    state: State,
    /// What NCBI's reader posted while it read the objects, not yet taken (`take_posts`).
    posts: Vec<Post>,
}

impl fmt::Debug for MaskFileReader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut debug = f.debug_struct("MaskFileReader");
        debug.field("path", &self.path);
        match &self.source {
            Source::Text(_) => debug.field("format", &MaskFileFormat::TextAsn),
            Source::Binary(rd) => debug
                .field("format", &MaskFileFormat::BinaryAsn)
                .field("offset", &rd.input().offset()),
        };
        debug.field("state", &self.state).finish()
    }
}

impl MaskFileReader {
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:933-947
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:160-181
    /// [`open_named`](Self::open_named) for a file that messages call by its path.
    #[cfg(test)]
    pub(crate) fn open(path: &Path) -> Result<Self, MaskFileError> {
        Self::open_named(path, path)
    }

    /// Opens the file and guesses its format from its first bytes; messages call the file
    /// `name` (the name on makeblastdb's command line). The text of the error for a file in
    /// neither format is NCBI's. A file that cannot be opened or read is a failed stream for the
    /// format guess (`x_TestInput`), which answers `eUnknown`: the same error. The file is then
    /// read in blocks as its objects are read (a read error ends the file, as the bad bit of
    /// NCBI's stream does).
    pub(crate) fn open_named(path: &Path, name: &Path) -> Result<Self, MaskFileError> {
        match std::fs::File::open(path) {
            Ok(file) => Self::from_reader(name, Box::new(file)),
            Err(_) => Err(unknown_encoding(name)),
        }
    }

    /// The reader of a file whose bytes are `data`; `path` is only used in messages.
    #[cfg(test)]
    pub(crate) fn from_bytes(path: &Path, data: Vec<u8>) -> Result<Self, MaskFileError> {
        Self::from_reader(path, Box::new(std::io::Cursor::new(data)))
    }

    /// The reader of the bytes `input` gives, read one byte per fill after the format guess (the
    /// bytes it keeps end at every byte); `path` is only used in messages.
    #[cfg(test)]
    pub(crate) fn from_read(path: &Path, input: Box<dyn Read>) -> Result<Self, MaskFileError> {
        let mut reader = Self::from_reader(path, input)?;
        match &mut reader.source {
            Source::Text(rd) => rd.input_mut().r.src.fill_one_byte_at_a_time(),
            Source::Binary(rd) => rd.input_mut().src.fill_one_byte_at_a_time(),
        }
        Ok(reader)
    }

    /// The reader of the file `input`; `path` is only used in messages.
    fn from_reader(path: &Path, mut input: Box<dyn Read>) -> Result<Self, MaskFileError> {
        let Some((format, sample)) = guess_format(&mut input) else {
            return Err(unknown_encoding(path));
        };
        let src = StreamSource::new(input, &sample);
        let source = match format {
            MaskFileFormat::TextAsn => {
                let stream = ObjectStream::text(src, Names::MASK);
                Source::Text(Box::new(stream))
            }
            MaskFileFormat::BinaryAsn => {
                let stream = ObjectStream::binary(src, Names::MASK);
                Source::Binary(Box::new(stream))
            }
        };
        Ok(MaskFileReader {
            path: path.to_path_buf(),
            source,
            state: State::Open,
            posts: Vec::new(),
        })
    }

    /// The bytes of the file that the reader keeps (the stream's buffer).
    #[cfg(test)]
    pub(crate) fn buffered(&self) -> usize {
        match &self.source {
            Source::Text(rd) => rd.input().buffered(),
            Source::Binary(rd) => rd.input().src.buffered(),
        }
    }

    #[cfg(test)]
    pub(crate) fn format(&self) -> MaskFileFormat {
        match self.source {
            Source::Text(_) => MaskFileFormat::TextAsn,
            Source::Binary(_) => MaskFileFormat::BinaryAsn,
        }
    }

    /// The next object of the file, as the root type `root`, read by `read` for either encoding;
    /// what was posted while it was read waits for `take_posts`.
    fn run<T>(
        &mut self,
        root: u16,
        text: impl FnOnce(&mut TypedReader<TextInput<FileSource>>) -> R<T>,
        binary: impl FnOnce(&mut TypedReader<BerReader<FileSource>>) -> R<T>,
    ) -> Result<T, ReadError> {
        let (result, mut posts) = match &mut self.source {
            Source::Text(stream) => stream.read_typed(root, text),
            Source::Binary(stream) => stream.read_typed(root, binary),
        };
        self.posts.append(&mut posts);
        result
    }

    /// The messages that NCBI's reader posted while it read the objects returned so far (and the
    /// one that failed), in the order posted: makeblastdb writes them as they come, before what
    /// it then does with the object or the error.
    pub(crate) fn take_posts(&mut self) -> Vec<Post> {
        std::mem::take(&mut self.posts)
    }

    /// A call after the reader stopped (makeblastdb does not read on after an error).
    fn stopped(&self) -> MaskFileError {
        if self.state == State::FailedEof {
            return MaskFileError::EofInContinuation {
                path: self.path.clone(),
            };
        }
        MaskFileError::Unsupported(format!(
            "mask data in '{}' after the error that stopped its reader",
            self.path.display()
        ))
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-263 (`CATCH_ALL`)
    /// The error of the reader's `fail`: what makeblastdb's `CATCH_ALL` reports for what NCBI's
    /// reader threw.
    fn fail(&mut self, fail: ReadError) -> MaskFileError {
        self.state = if matches!(fail, ReadError::Eof) {
            State::FailedEof
        } else {
            State::Failed
        };
        match fail {
            ReadError::Serial(e) | ReadError::Exception(e) => {
                MaskFileError::Ncbi(WriteDbError::Exception(e))
            }
            ReadError::OutOfMemory => MaskFileError::Ncbi(WriteDbError::OutOfMemory),
            ReadError::Std(what) => MaskFileError::Ncbi(WriteDbError::Other(what.to_string())),
            ReadError::Unsupported(what) => MaskFileError::Unsupported(what),
            ReadError::Eof => MaskFileError::EofInContinuation {
                path: self.path.clone(),
            },
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:947-960
    /// The next `Blast-db-mask-info`; `None` at the end of the file. A `CEofException` of the
    /// serial reader (nothing left but white space and comments, or a file that ends inside the
    /// type name or the first bytes of an object) is the end, as in `x_ProcessMaskData`; a
    /// partly read object is dropped.
    pub(crate) fn next_info(&mut self) -> Result<Option<MaskInfo>, MaskFileError> {
        match self.state {
            State::Open => {}
            State::Ended => return Ok(None),
            State::Failed | State::FailedEof => return Err(self.stopped()),
        }
        match self.run(named::BLAST_DB_MASK_INFO, Build::info, Build::info) {
            Ok(info) => Ok(Some(info)),
            Err(ReadError::Eof) => {
                self.state = State::Ended;
                Ok(None)
            }
            Err(fail) => Err(self.fail(fail)),
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:994-998
    /// The `Blast-mask-list` that follows an object with `more TRUE`. A file that ends there is
    /// [`MaskFileError::EofInContinuation`]. After an error every call is an error again.
    pub(crate) fn next_list(&mut self) -> Result<MaskList, MaskFileError> {
        match self.state {
            State::Open => {}
            State::Ended => {
                return Err(MaskFileError::EofInContinuation {
                    path: self.path.clone(),
                })
            }
            State::Failed | State::FailedEof => return Err(self.stopped()),
        }
        match self.run(named::BLAST_MASK_LIST, Build::list, Build::list) {
            Ok(list) => Ok(list),
            Err(fail) => Err(self.fail(fail)),
        }
    }
}
