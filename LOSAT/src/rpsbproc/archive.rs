//! The reader of rpsbproc's input: ASN.1 text `Blast4-archive` objects, one after another
//! from one stream (rpsblast `-outfmt 11` writes one per query batch), read one object at a time
//! as NCBI's rpsbproc reads them.
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:780-792
//!
//! An object is read whole before rpsbproc uses it, and any failure of NCBI's typed reader (a
//! member or variant that the module does not have, a member out of order or twice, a missing
//! mandatory member, a `,` before `}`, a value of the wrong form, the end of the text inside the
//! object, anything after the last object that is not an object) ends the input: the objects
//! before it are processed and the run ends normally. The grammar below follows the ASN.1
//! modules for every type that it descends into (names, order, OPTIONAL and DEFAULT members,
//! CHOICE variants, ENUMERATED values); the values that rpsbproc does not use are read by the
//! typed reader of the type table (`ArchiveStream::typed_member`: the request bodies other than
//! `queue-search`, the `Seqdesc` variants other than `title` and `user`, `Seq-annot`, `Int-fuzz`,
//! the `Seq-loc` variants that the grammar does not read, the `Seq-align` members `bounds` and
//! `ext`, the segments other than `denseg`, `std` and `packed`, the large `Blast4-value`
//! variants, the `phi-alignments`, `pssm` and `simple-results` of the results, the `Seq-ext`
//! variants other than `delta` and the `loc` items of a `delta`, `Seq-hist` and `gap` data), over
//! the same lexical reader.
//!
//! A byte outside 0x20..0x7E in a VisibleString is read as `#`, with NCBI's post of it
//! (`Bad char`, 803.7, the stack trace of the string), as NCBI's reader does by default
//! (`eFNP_ReplaceAndWarn` of the parameter `[SERIAL] WRONG_CHARS_READ`, a registry entry, which
//! the application start rejects in a file and in `NCBI_CONFIG__SERIAL__WRONG_CHARS_READ`); a
//! UTF8String keeps its bytes. The grammar passes the frames of NCBI's object stack down to the
//! reader of each value (`asn_text_reader::Frame`: the member, variant or element that holds it),
//! which give the trace of its strings and the frames below which the typed reader reads. The
//! posts of an object come before rpsbproc uses it (`ArchiveRead`).
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1897
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:291-308

use std::io::Read;

use crate::report::serial::asn_input::{Input, Names, TextInput};
use crate::report::serial::asn_text_reader::{
    AsnTextError, AsnTextErrorKind, Frame, Reader, Result, StreamSource,
};
use crate::report::serial::asn_types::read::{ObjectStream, Post, ReadError};
use crate::report::serial::asn_types::{self, table::named};
use crate::report::serial::pssm::{SEQ_INST_MOL, SEQ_INST_REPR};
use crate::seqid::model::{ObjectId, ParsedSeqId};

/// One `Blast4-archive`, with what rpsbproc takes from it.
// NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:67-74
#[derive(Debug)]
pub(crate) struct Archive {
    /// `request.ident`.
    pub(crate) ident: Option<Vec<u8>>,
    /// `request.body`.
    pub(crate) body: RequestBody,
    pub(crate) results: Results,
}

/// `Blast4-request-body`: `queue-search`, or another variant.
#[derive(Debug)]
pub(crate) enum RequestBody {
    QueueSearch(QueueSearch),
    Other,
}

/// `Blast4-queue-search-request`.
#[derive(Debug)]
pub(crate) struct QueueSearch {
    pub(crate) program: Vec<u8>,
    pub(crate) queries: Queries,
    pub(crate) subject: Subject,
    /// Whether `algorithm-options` is set.
    pub(crate) has_algorithm_options: bool,
    /// The value of the first `algorithm-options` parameter named `MatrixName`
    /// (`CBlast4_parameters::GetParamByName`).
    pub(crate) matrix_name: Option<ParamValue>,
    /// The value of the first parameter named `EvalueThreshold`.
    pub(crate) evalue_threshold: Option<ParamValue>,
}

/// `Blast4-queries`: the query Bioseqs of a `bioseq-set`, or another variant.
#[derive(Debug)]
pub(crate) enum Queries {
    BioseqSet(Vec<SeqEntry>),
    Other,
}

/// `Blast4-subject`: a database name, or the name of another variant.
#[derive(Debug)]
pub(crate) enum Subject {
    Database(Vec<u8>),
    Other(&'static str),
}

/// The value of a `Blast4-parameter` as far as rpsbproc reads it.
#[derive(Debug)]
pub(crate) enum ParamValue {
    String(Vec<u8>),
    Cutoff(Cutoff),
    /// Another `Blast4-value` variant, by name.
    Other(&'static str),
}

/// `Blast4-cutoff`.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Cutoff {
    EValue(f64),
    RawScore,
}

/// A `Seq-entry` of the query set: a Bioseq, or a nested `set`.
#[derive(Debug)]
pub(crate) enum SeqEntry {
    Seq(Query),
    Set,
}

/// A query `Bioseq`, with what `FillSequenceFromBioseq` takes from it.
#[derive(Debug, Default)]
pub(crate) struct Query {
    pub(crate) ids: Vec<ParsedSeqId>,
    /// The first `title` descriptor.
    pub(crate) title: Option<Vec<u8>>,
    /// The first descriptor that LOSAT does not read for this query: a `CFastaReader` user
    /// object whose `DefLine` field is not a string (`ObtainUserDeflineFromBioseq` asks for its
    /// string).
    pub(crate) unread_descriptor: Option<&'static str>,
    /// Whether `annot` is set (`GetDescrTitleFromBioseq` looks there for a title when the
    /// descriptors have none).
    pub(crate) has_annot: bool,
    /// `inst.mol`.
    pub(crate) mol: i32,
    /// `inst.length`.
    pub(crate) length: Option<u32>,
    /// `inst.seq-data`: its variant and its letters (or bytes).
    pub(crate) seq_data: Option<(&'static str, Vec<u8>)>,
    /// The `seq-data` of the literals of `inst.ext delta`, in order (literals without data or
    /// with `gap` data and `loc` items are left out, as `Get1LetterSeqData` skips them).
    pub(crate) delta_literals: Vec<(&'static str, Vec<u8>)>,
}

/// `Blast4-get-search-results-reply`, with what rpsbproc reads.
#[derive(Debug, Default)]
pub(crate) struct Results {
    pub(crate) alignments: Option<Vec<Align>>,
    pub(crate) masks: Option<Vec<Mask>>,
}

/// A `Seq-align` of the results.
#[derive(Debug)]
pub(crate) struct Align {
    pub(crate) segs: Segs,
    /// `ParseAlignScores`: the values of the scores named `score`, `e_value`, `bit_score` and
    /// `num_ident` (the last one of each name).
    pub(crate) score: Option<i32>,
    pub(crate) e_value: Option<f64>,
    pub(crate) bit_score: Option<f64>,
    pub(crate) num_ident: Option<i32>,
    /// The first of those scores whose value is not of the type rpsbproc asks for
    /// (`CScore::C_Value::GetInt` or `GetReal` would throw).
    pub(crate) score_of_other_type: Option<&'static str>,
}

/// `Seq-align.segs`: a `Dense-seg`, the `Std-seg`s of a translated search (rpstblastn), or the
/// name of another variant.
#[derive(Debug)]
pub(crate) enum Segs {
    Denseg(DenseSeg),
    /// (Boxed: its two Seq-ids would make every `Segs` four times as large.)
    Std(Box<StdAlign>),
    Other(&'static str),
}

/// The `Std-seg`s of a Seq-align as far as rpsbproc reads them: the query row's id that
/// `CSeq_align::GetSeq_id` finds for index 0, the second id of the first segment
/// (`ParseAlignSegs`'s hit id) and the first two locations of each segment. The ids of the other
/// places are read and left (an archive of a genome holds hundreds of thousands of segments).
// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:94-98
#[derive(Debug)]
pub(crate) struct StdAlign {
    pub(crate) query_id: QueryRowId,
    /// `ids[1]` of the first segment; `None` when it has fewer than 2 ids.
    pub(crate) hit_id: Option<ParsedSeqId>,
    pub(crate) segs: Vec<StdPair>,
}

/// What `CSeq_align::GetSeq_id` finds for index 0 in the segments, in order: the first id of the first
/// segment with ids, else the id of a first location of a single id (`null` has none and goes
/// on to the next segment).
#[derive(Debug)]
pub(crate) enum QueryRowId {
    Id(ParsedSeqId),
    /// A first location of another variant (`packed-int`, `mix`, ...), by name, whose id LOSAT
    /// does not look up.
    Location(&'static str),
    /// No segment gives an id.
    Missing,
}

/// The first two locations of a `Std-seg` (the query's, the subject's).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct StdPair {
    pub(crate) query: StdLoc,
    pub(crate) hit: StdLoc,
}

/// A `Seq-loc` of a `Std-seg`, as far as `ParseAlignSegs` reads it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum StdLoc {
    /// The segment has no location at this place (fewer than 2).
    Missing,
    Null,
    Empty,
    Whole,
    /// A `Seq-interval`; `from` and `to` are `TSeqPos` (unsigned), `strand` the value of its
    /// `Na-strand` when it is set.
    Int {
        from: u32,
        to: u32,
        strand: Option<i32>,
    },
    /// A `Seq-point`.
    Pnt {
        point: u32,
        strand: Option<i32>,
    },
    /// Another variant (`packed-int`, `packed-pnt`, `mix`, `equiv`, `bond`, `feat`).
    Other,
}

/// A `Seq-interval` as read.
struct Interval {
    from: u32,
    to: u32,
    strand: Option<i32>,
    id: ParsedSeqId,
}

/// `Dense-seg`.
#[derive(Debug)]
pub(crate) struct DenseSeg {
    pub(crate) ids: Vec<ParsedSeqId>,
    pub(crate) starts: Vec<i32>,
    pub(crate) lens: Vec<u32>,
}

/// The first location of a `Blast4-mask` (the one `ConvertBlastMaskListToSeqLocInfoVector`
/// reads).
#[derive(Debug)]
pub(crate) enum MaskLocation {
    /// A `packed-int`: the ids of its intervals.
    PackedInt(Vec<ParsedSeqId>),
    /// Another `Seq-loc` variant, by name.
    Other(&'static str),
}

/// A `Blast4-mask` with locations (a mask without any is not used).
#[derive(Debug)]
pub(crate) struct Mask {
    pub(crate) first: Option<MaskLocation>,
}

/// The objects of one input stream: the hand grammar below over the lexical reader of the typed
/// reader's input (`ObjectStream::input_mut`), which reads the values that the grammar does not
/// take (`typed_member`).
pub(crate) struct ArchiveStream<R: Read> {
    objects: ObjectStream<TextInput<StreamSource<R>>>,
    /// What NCBI's reader posts while it reads the object, in order: the hand grammar's posts
    /// (`Trail::replaced`) and the typed reader's.
    posts: Vec<Post>,
    done: bool,
}

/// What reading one object gives (`ArchiveStream::next_archive`): the object, the end of the input
/// (`Ok(None)`) or LOSAT's refusal, and the messages that NCBI's reader posts while it reads the
/// object, which come before rpsbproc uses it (also when the object then fails to read).
pub(crate) type ArchiveRead = (
    std::result::Result<Option<Archive>, AsnTextError>,
    Vec<Post>,
);

/// The name of the root type (`SkipFileHeader`'s check, the first name of the stack traces).
fn blast4_archive() -> &'static str {
    asn_types::type_info(named::BLAST4_ARCHIVE).name()
}

/// The `Blast4-request-body` variants.
// NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:76-98
const REQUEST_BODY: [&str; 15] = [
    "finish-params",
    "get-databases",
    "get-matrices",
    "get-parameters",
    "get-paramsets",
    "get-programs",
    "get-search-results",
    "get-sequences",
    "queue-search",
    "get-request-info",
    "get-sequence-parts",
    "get-windowmasked-taxids",
    "get-protocol-info",
    "get-search-info",
    "get-databases-ex",
];

/// The `Blast4-value` variants.
// NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:512-551
const BLAST4_VALUE: [&str; 27] = [
    "big-integer",
    "bioseq",
    "boolean",
    "cutoff",
    "integer",
    "matrix",
    "real",
    "seq-align",
    "seq-id",
    "seq-loc",
    "strand-type",
    "string",
    "big-integer-list",
    "bioseq-list",
    "boolean-list",
    "cutoff-list",
    "integer-list",
    "matrix-list",
    "real-list",
    "seq-align-list",
    "seq-id-list",
    "seq-loc-list",
    "strand-type-list",
    "string-list",
    "bioseq-set",
    "seq-align-set",
    "query-mask",
];

// NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:493-497
const BLAST4_STRAND_TYPE: [(&str, i32); 3] = [
    ("forward-strand", 1),
    ("reverse-strand", 2),
    ("both-strands", 3),
];

// NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:423-431
const BLAST4_FRAME_TYPE: [(&str, i32); 7] = [
    ("notset", 0),
    ("plus1", 1),
    ("plus2", 2),
    ("plus3", 3),
    ("minus1", 4),
    ("minus2", 5),
    ("minus3", 6),
];

// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:115-121
const NA_STRAND: [(&str, i32); 6] = [
    ("unknown", 0),
    ("plus", 1),
    ("minus", 2),
    ("both", 3),
    ("both-rev", 4),
    ("other", 255),
];

// NCBI reference (598d8ae6): c++/src/objects/seqset/seqset.asn:26-51
const BIOSEQ_SET_CLASS: [(&str, i32); 25] = [
    ("not-set", 0),
    ("nuc-prot", 1),
    ("segset", 2),
    ("conset", 3),
    ("parts", 4),
    ("gibb", 5),
    ("gi", 6),
    ("genbank", 7),
    ("pir", 8),
    ("pub-set", 9),
    ("equiv", 10),
    ("swissprot", 11),
    ("pdb-entry", 12),
    ("mut-set", 13),
    ("pop-set", 14),
    ("phy-set", 15),
    ("eco-set", 16),
    ("gen-prod-set", 17),
    ("wgs-set", 18),
    ("named-annot", 19),
    ("named-annot-prod", 20),
    ("read-set", 21),
    ("paired-end-reads", 22),
    ("small-genome-set", 23),
    ("other", 255),
];

// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:283-294 (Seq-inst topology, strand)
const SEQ_INST_TOPOLOGY: [(&str, i32); 5] = [
    ("not-set", 0),
    ("linear", 1),
    ("circular", 2),
    ("tandem", 3),
    ("other", 255),
];
const SEQ_INST_STRAND: [(&str, i32); 5] = [
    ("not-set", 0),
    ("ss", 1),
    ("ds", 2),
    ("mixed", 3),
    ("other", 255),
];

// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:47-73
const SEQDESC: [&str; 25] = [
    "mol-type",
    "modif",
    "method",
    "name",
    "title",
    "org",
    "comment",
    "num",
    "maploc",
    "pir",
    "genbank",
    "pub",
    "region",
    "user",
    "sp",
    "dbxref",
    "embl",
    "create-date",
    "update-date",
    "prf",
    "pdb",
    "het",
    "source",
    "molinfo",
    "modelev",
];

// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:346-358
// The letter codes (`IUPACna`, `IUPACaa`, `NCBIeaa`) are `StringStore` and the others
// `OCTET STRING` (seq.asn:400-461).
const SEQ_DATA: [&str; 11] = [
    "iupacna",
    "iupacaa",
    "ncbi2na",
    "ncbi4na",
    "ncbi8na",
    "ncbipna",
    "ncbi8aa",
    "ncbieaa",
    "ncbipaa",
    "ncbistdaa",
    "gap",
];

// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:79-90
const SEQ_LOC: [&str; 11] = [
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
];

// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:32-40
const SEQ_ALIGN_SEGS: [&str; 7] = [
    "dendiag", "denseg", "std", "packed", "disc", "spliced", "sparse",
];

/// The type of `Seq-align.segs` (a CHOICE without a name of its own) in the type table.
fn seq_align_segs() -> u16 {
    let members = asn_types::members(named::SEQ_ALIGN);
    members[asn_types::member_index(named::SEQ_ALIGN, "segs")].ty
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:22-29
const SEQ_ALIGN_TYPE: [(&str, i32); 6] = [
    ("not-set", 0),
    ("global", 1),
    ("diags", 2),
    ("partial", 3),
    ("disc", 4),
    ("other", 255),
];

impl<R: Read> ArchiveStream<R> {
    /// The objects of `input`, whose first bytes `first` were already taken from it (the bytes
    /// after the white space that rpsbproc skips to tell ASN.1 from XML).
    pub(crate) fn new(input: R, first: &[u8]) -> Self {
        let mut objects = ObjectStream::text(StreamSource::new(input, first), Names::ARCHIVE);
        objects.input_mut().r.trail.root = blast4_archive();
        ArchiveStream {
            objects,
            posts: Vec::new(),
            done: false,
        }
    }

    /// The lexical reader.
    #[inline]
    fn r(&mut self) -> &mut Reader<StreamSource<R>> {
        &mut self.objects.input_mut().r
    }

    /// The bytes the stream keeps.
    #[cfg(test)]
    pub(crate) fn buffered(&self) -> usize {
        self.objects.input().r.src.buffered()
    }

    /// Reads one byte per fill from now on (`StreamSource::fill_one_byte_at_a_time`).
    #[cfg(test)]
    pub(crate) fn fill_one_byte_at_a_time(&mut self) {
        self.r().src.fill_one_byte_at_a_time();
    }

    /// The next object, with the messages posted while it was read: `Ok(None)` when NCBI's
    /// reader fails at it (the end of the input), the object when it reads, or the refusal of
    /// what LOSAT does not read in it.
    pub(crate) fn next_archive(&mut self) -> ArchiveRead {
        if self.done {
            return (Ok(None), Vec::new());
        }
        let result = match self.archive_object() {
            Ok(archive) => Ok(Some(archive)),
            Err(error) => {
                self.done = true;
                match error.kind {
                    AsnTextErrorKind::Format | AsnTextErrorKind::Eof => Ok(None),
                    AsnTextErrorKind::Unsupported => Err(error),
                }
            }
        };
        self.take_replaced();
        (result, std::mem::take(&mut self.posts))
    }

    /// The posts of the hand grammar's strings so far, after those posted before them.
    fn take_replaced(&mut self) {
        let replaced = std::mem::take(&mut self.r().trail.replaced);
        self.posts.extend(replaced.into_iter().map(Post::BadChar));
    }

    /// The value of the member or variant of the frame `at`, of the SEQUENCE or CHOICE `parent` of
    /// the type table, read by the typed reader below the frames of `at` as NCBI's reader reads
    /// it: rpsbproc does not use it, and NCBI's reader fails, posts and reads on in it as in any
    /// value of its object. Its failure ends the input (NCBI's rpsbproc stops at any exception of
    /// its reader), except LOSAT's refusal.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:780-792
    #[cold]
    #[inline(never)]
    fn typed_member(&mut self, parent: u16, at: &Frame<'_>) -> Result<()> {
        self.take_replaced();
        let name = at.name().expect("a member or variant is being read");
        let ty = asn_types::members(parent)[asn_types::member_index(parent, name)].ty;
        let root = self.r().trail.root;
        let (value, posts) = self.objects.read_value(ty, root, &at.steps());
        self.posts.extend(posts);
        match value {
            Ok(_) => Ok(()),
            Err(ReadError::Unsupported(what)) => Err(self.r().not_supported(&format!(
                "an archive value that LOSAT's reader does not read ({what})"
            ))),
            Err(_) => Err(self.r().error("the typed reader failed")),
        }
    }

    fn real(&mut self) -> Result<f64> {
        self.r().real()
    }

    fn string(&mut self, at: &Frame<'_>) -> Result<Vec<u8>> {
        self.r().string(at)
    }

    /// `ReadFileHeader` and the type check of the root reader, then the object. Each object is
    /// read by a new object stream of NCBI's (`ObjStreamIn`), whose line count starts again at the
    /// byte after the last one the object before it took (`ObjectStream`, `Input::new_object`).
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:892-898
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/objutils.hpp:75-83
    fn archive_object(&mut self) -> Result<Archive> {
        self.objects.input_mut().new_object();
        let r = self.r();
        let id = r.file_header()?;
        if !id.is_empty() && r.bytes(id) != r.trail.root.as_bytes() {
            return Err(r.error("incompatible type"));
        }
        let mut members = r.begin_members(&Frame::ROOT, &["request", "results", "messages"])?;
        let (mut request, mut results) = (None, None);
        while let Some(index) = self.r().next_member(&mut members)? {
            match index {
                0 => request = Some(self.request(members.at())?),
                1 => results = Some(self.results(members.at())?),
                _ => self.messages(members.at())?,
            }
        }
        let (Some((ident, body)), Some(results)) = (request, results) else {
            return Err(self.r().error("a mandatory member is missing"));
        };
        Ok(Archive {
            ident,
            body,
            results,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:56-63
    fn request(&mut self, at: &Frame<'_>) -> Result<(Option<Vec<u8>>, RequestBody)> {
        let mut members = self.r().begin_members(at, &["ident", "body"])?;
        let (mut ident, mut body) = (None, None);
        while let Some(index) = self.r().next_member(&mut members)? {
            match index {
                0 => ident = Some(self.string(members.at())?),
                _ => {
                    let variant = self.r().choice(members.at(), &REQUEST_BODY)?;
                    body = Some(if variant.index() == 8 {
                        RequestBody::QueueSearch(self.queue_search(&variant)?)
                    } else {
                        self.typed_member(named::BLAST4_REQUEST_BODY, &variant)?;
                        RequestBody::Other
                    });
                }
            }
        }
        let body = body.ok_or_else(|| self.r().error("member body missing"))?;
        Ok((ident, body))
    }

    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:156-166
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:138-142
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:499-503
    fn queue_search(&mut self, at: &Frame<'_>) -> Result<QueueSearch> {
        const NAMES: [&str; 8] = [
            "program",
            "service",
            "queries",
            "subject",
            "paramset",
            "algorithm-options",
            "program-options",
            "format-options",
        ];
        let mut members = self.r().begin_members(at, &NAMES)?;
        let (mut program, mut service, mut queries, mut subject) = (None, None, None, None);
        let mut search = QueueSearch {
            program: Vec::new(),
            queries: Queries::Other,
            subject: Subject::Other("sequences"),
            has_algorithm_options: false,
            matrix_name: None,
            evalue_threshold: None,
        };
        while let Some(index) = self.r().next_member(&mut members)? {
            let at = members.at();
            match index {
                0 => program = Some(self.string(at)?),
                1 => service = Some(self.string(at)?),
                2 => {
                    const QUERIES: [&str; 3] = ["pssm", "seq-loc-list", "bioseq-set"];
                    let variant = self.r().choice(at, &QUERIES)?;
                    queries = Some(if variant.index() == 2 {
                        Queries::BioseqSet(self.bioseq_set(&variant)?)
                    } else {
                        self.typed_member(named::BLAST4_QUERIES, &variant)?;
                        Queries::Other
                    });
                }
                3 => {
                    const SUBJECT: [&str; 3] = ["database", "sequences", "seq-loc-list"];
                    let variant = self.r().choice(at, &SUBJECT)?;
                    subject = Some(if variant.index() == 0 {
                        Subject::Database(self.string(&variant)?)
                    } else {
                        self.typed_member(named::BLAST4_SUBJECT, &variant)?;
                        Subject::Other(SUBJECT[variant.index()])
                    });
                }
                4 => {
                    self.string(at)?;
                }
                5 => {
                    search.has_algorithm_options = true;
                    let (matrix, evalue) = self.parameters(at)?;
                    search.matrix_name = matrix;
                    search.evalue_threshold = evalue;
                }
                _ => {
                    self.parameters(at)?;
                }
            }
        }
        let (Some(program), Some(_), Some(queries), Some(subject)) =
            (program, service, queries, subject)
        else {
            return Err(self.r().error("a mandatory member is missing"));
        };
        search.program = program;
        search.queries = queries;
        search.subject = subject;
        Ok(search)
    }

    /// `Blast4-parameters`: the values of the first parameters named `MatrixName` and
    /// `EvalueThreshold`.
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:464-467,505
    // NCBI reference (598d8ae6): c++/src/objects/blast/Blast4_parameters.cpp:186-198
    fn parameters(&mut self, at: &Frame<'_>) -> Result<(Option<ParamValue>, Option<ParamValue>)> {
        let (mut matrix, mut evalue) = (None, None);
        let mut items = self.r().begin_list(at)?;
        while self.r().next_item(&mut items)? {
            let mut members = self.r().begin_members(items.at(), &["name", "value"])?;
            let (mut name, mut value) = (None, None);
            while let Some(index) = self.r().next_member(&mut members)? {
                match index {
                    0 => name = Some(self.string(members.at())?),
                    _ => value = Some(self.blast4_value(members.at())?),
                }
            }
            let (Some(name), Some(value)) = (name, value) else {
                return Err(self.r().error("a mandatory member is missing"));
            };
            if name == b"MatrixName" && matrix.is_none() {
                matrix = Some(value);
            } else if name == b"EvalueThreshold" && evalue.is_none() {
                evalue = Some(value);
            }
        }
        self.r().end_block()?;
        Ok((matrix, evalue))
    }

    /// A `Blast4-value` (the variants listed at `BLAST4_VALUE`).
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:375-378
    fn blast4_value(&mut self, at: &Frame<'_>) -> Result<ParamValue> {
        let variant = self.r().choice(at, &BLAST4_VALUE)?;
        let name = BLAST4_VALUE[variant.index()];
        match name {
            "big-integer" => {
                self.r().int8()?;
            }
            "boolean" => {
                self.r().boolean()?;
            }
            "cutoff" => {
                return Ok(ParamValue::Cutoff(
                    if self
                        .r()
                        .choice(&variant, &["e-value", "raw-score"])?
                        .index()
                        == 0
                    {
                        Cutoff::EValue(self.real()?)
                    } else {
                        self.r().int4()?;
                        Cutoff::RawScore
                    },
                ));
            }
            "integer" => {
                self.r().int4()?;
            }
            "real" => {
                self.real()?;
            }
            "strand-type" => {
                self.r().enumerated(&BLAST4_STRAND_TYPE)?;
            }
            "string" => return Ok(ParamValue::String(self.string(&variant)?)),
            _ => self.typed_member(named::BLAST4_VALUE, &variant)?,
        }
        Ok(ParamValue::Other(name))
    }

    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:287-299
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:433-438,454-457
    fn results(&mut self, at: &Frame<'_>) -> Result<Results> {
        const NAMES: [&str; 7] = [
            "alignments",
            "phi-alignments",
            "masks",
            "ka-blocks",
            "search-stats",
            "pssm",
            "simple-results",
        ];
        let mut results = Results::default();
        let mut members = self.r().begin_members(at, &NAMES)?;
        while let Some(index) = self.r().next_member(&mut members)? {
            let at = members.at();
            match index {
                0 => results.alignments = Some(self.seq_align_set(at)?),
                2 => {
                    let mut masks = Vec::new();
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        masks.push(self.blast4_mask(items.at())?);
                    }
                    self.r().end_block()?;
                    results.masks = Some(masks);
                }
                3 => {
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        let mut members = self
                            .r()
                            .begin_members(items.at(), &["lambda", "k", "h", "gapped"])?;
                        let mut seen = 0;
                        while let Some(index) = self.r().next_member(&mut members)? {
                            seen += 1;
                            if index == 3 {
                                self.r().boolean()?;
                            } else {
                                self.real()?;
                            }
                        }
                        if seen != 4 {
                            return Err(self.r().error("a mandatory member is missing"));
                        }
                    }
                    self.r().end_block()?;
                }
                4 => {
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        self.string(items.at())?;
                    }
                    self.r().end_block()?;
                }
                _ => self.typed_member(named::BLAST4_GET_SEARCH_RESULTS_REPLY, at)?,
            }
        }
        Ok(results)
    }

    fn blast4_mask(&mut self, at: &Frame<'_>) -> Result<Mask> {
        let mut members = self.r().begin_members(at, &["locations", "frame"])?;
        let (mut first, mut locations, mut frame) = (None, false, false);
        while let Some(index) = self.r().next_member(&mut members)? {
            if index == 0 {
                locations = true;
                let mut items = self.r().begin_list(members.at())?;
                while self.r().next_item(&mut items)? {
                    let location = self.mask_location(items.at())?;
                    if first.is_none() {
                        first = Some(location);
                    }
                }
                self.r().end_block()?;
            } else {
                frame = true;
                self.r().enumerated(&BLAST4_FRAME_TYPE)?;
            }
        }
        if !(locations && frame) {
            return Err(self.r().error("a mandatory member is missing"));
        }
        Ok(Mask { first })
    }

    /// A `Seq-loc` of a mask: the ids of a `packed-int`, the name of another variant.
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:93-101
    // `from` and `to` are `TSeqPos` (unsigned) in NCBI's classes.
    fn mask_location(&mut self, at: &Frame<'_>) -> Result<MaskLocation> {
        let variant = self.r().choice(at, &SEQ_LOC)?;
        match SEQ_LOC[variant.index()] {
            "int" => {
                self.seq_interval(&variant)?;
                Ok(MaskLocation::Other("int"))
            }
            "packed-int" => {
                let mut ids = Vec::new();
                let mut items = self.r().begin_list(&variant)?;
                while self.r().next_item(&mut items)? {
                    ids.push(self.seq_interval(items.at())?);
                }
                self.r().end_block()?;
                Ok(MaskLocation::PackedInt(ids))
            }
            name => {
                self.typed_member(named::SEQ_LOC, &variant)?;
                Ok(MaskLocation::Other(name))
            }
        }
    }

    /// A `Seq-interval`: its id.
    fn seq_interval(&mut self, at: &Frame<'_>) -> Result<ParsedSeqId> {
        Ok(self.interval(at)?.id)
    }

    /// A `Seq-interval`.
    fn interval(&mut self, at: &Frame<'_>) -> Result<Interval> {
        const NAMES: [&str; 6] = ["from", "to", "strand", "id", "fuzz-from", "fuzz-to"];
        let mut members = self.r().begin_members(at, &NAMES)?;
        let (mut from, mut to, mut strand, mut id) = (None, None, None, None);
        while let Some(index) = self.r().next_member(&mut members)? {
            match index {
                0 => from = Some(self.r().uint4()?),
                1 => to = Some(self.r().uint4()?),
                2 => strand = Some(self.r().enumerated(&NA_STRAND)?.1),
                3 => id = Some(self.r().seq_id(members.at())?),
                _ => self.typed_member(named::SEQ_INTERVAL, members.at())?,
            }
        }
        match (from, to, id) {
            (Some(from), Some(to), Some(id)) => Ok(Interval {
                from,
                to,
                strand,
                id,
            }),
            _ => Err(self.r().error("a mandatory member is missing")),
        }
    }

    /// A `Packed-seg`, read for its grammar (`CSeq_align::GetSeq_id` throws for it before
    /// rpsbproc reads more of it, `query::query_row_id`).
    // NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:83-92
    // `starts` and `lens` are `TSeqPos` (unsigned) in NCBI's class (seqalign.def:28-30).
    fn packed_seg(&mut self, at: &Frame<'_>) -> Result<()> {
        const NAMES: [&str; 8] = [
            "dim", "numseg", "ids", "starts", "present", "lens", "strands", "scores",
        ];
        let mut seen = [false; 8];
        let mut members = self.r().begin_members(at, &NAMES)?;
        while let Some(index) = self.r().next_member(&mut members)? {
            seen[index] = true;
            match index {
                0 | 1 => {
                    self.r().int4()?;
                }
                4 => {
                    self.r().octets()?;
                }
                _ => {
                    let mut items = self.r().begin_list(members.at())?;
                    while self.r().next_item(&mut items)? {
                        match index {
                            2 => {
                                self.r().seq_id(items.at())?;
                            }
                            3 | 5 => {
                                self.r().uint4()?;
                            }
                            6 => {
                                self.r().enumerated(&NA_STRAND)?;
                            }
                            _ => {
                                let mut scratch = Align {
                                    segs: Segs::Other("dendiag"),
                                    score: None,
                                    e_value: None,
                                    bit_score: None,
                                    num_ident: None,
                                    score_of_other_type: None,
                                };
                                self.align_score(items.at(), &mut scratch)?;
                            }
                        }
                    }
                    self.r().end_block()?;
                }
            }
        }
        if !(seen[1] && seen[2] && seen[3] && seen[4] && seen[5]) {
            return Err(self.r().error("a mandatory member is missing"));
        }
        Ok(())
    }

    /// The `std` segments of a Seq-align (SEQUENCE OF Std-seg).
    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:339-358
    // (`ids` comes before `loc` in the module's order, so a location's id is kept only for the
    // first location of a segment without ids while no id is found.)
    fn std_segs(&mut self, at: &Frame<'_>) -> Result<StdAlign> {
        const NAMES: [&str; 4] = ["dim", "ids", "loc", "scores"];
        let mut align = StdAlign {
            query_id: QueryRowId::Missing,
            hit_id: None,
            segs: Vec::new(),
        };
        let mut items = self.r().begin_list(at)?;
        while self.r().next_item(&mut items)? {
            let first = align.segs.is_empty();
            let mut pair = StdPair {
                query: StdLoc::Missing,
                hit: StdLoc::Missing,
            };
            let mut ids_found = false;
            let mut loc = false;
            let mut members = self.r().begin_members(items.at(), &NAMES)?;
            while let Some(index) = self.r().next_member(&mut members)? {
                match index {
                    0 => {
                        self.r().int4()?;
                    }
                    1 => {
                        let mut id_items = self.r().begin_list(members.at())?;
                        let mut count = 0usize;
                        while self.r().next_item(&mut id_items)? {
                            let id = self.r().seq_id(id_items.at())?;
                            if count == 0 && matches!(align.query_id, QueryRowId::Missing) {
                                align.query_id = QueryRowId::Id(id);
                            } else if count == 1 && first {
                                align.hit_id = Some(id);
                            }
                            count += 1;
                        }
                        ids_found = count > 0;
                        self.r().end_block()?;
                    }
                    2 => {
                        loc = true;
                        let mut loc_items = self.r().begin_list(members.at())?;
                        let mut count = 0usize;
                        while self.r().next_item(&mut loc_items)? {
                            let wanted = count == 0
                                && !ids_found
                                && matches!(align.query_id, QueryRowId::Missing);
                            let (location, id, name) = self.std_loc(loc_items.at(), wanted)?;
                            if wanted {
                                align.query_id = match (location, id) {
                                    (_, Some(id)) => QueryRowId::Id(id),
                                    (StdLoc::Null, None) => QueryRowId::Missing,
                                    (_, None) => QueryRowId::Location(name),
                                };
                            }
                            match count {
                                0 => pair.query = location,
                                1 => pair.hit = location,
                                _ => {}
                            }
                            count += 1;
                        }
                        self.r().end_block()?;
                    }
                    _ => {
                        let mut scratch = Align {
                            segs: Segs::Other("dendiag"),
                            score: None,
                            e_value: None,
                            bit_score: None,
                            num_ident: None,
                            score_of_other_type: None,
                        };
                        let mut score_items = self.r().begin_list(members.at())?;
                        while self.r().next_item(&mut score_items)? {
                            self.align_score(score_items.at(), &mut scratch)?;
                        }
                        self.r().end_block()?;
                    }
                }
            }
            if !loc {
                return Err(self.r().error("a mandatory member is missing"));
            }
            align.segs.push(pair);
        }
        self.r().end_block()?;
        Ok(align)
    }

    /// A `Seq-loc` of a `Std-seg` and its variant's name: `null`, `empty`, `whole`, `int` and
    /// `pnt` read, the other variants read as values; with `keep_id`, the id of a location of a
    /// single id (`empty`, `whole`, `int`, `pnt`).
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:103-107
    // `point` is `TSeqPos` (unsigned) in NCBI's class.
    fn std_loc(
        &mut self,
        at: &Frame<'_>,
        keep_id: bool,
    ) -> Result<(StdLoc, Option<ParsedSeqId>, &'static str)> {
        let variant = self.r().choice(at, &SEQ_LOC)?;
        let name = SEQ_LOC[variant.index()];
        let keep = |id: ParsedSeqId| keep_id.then_some(id);
        let (location, id) = match name {
            "null" => {
                self.r().null()?;
                (StdLoc::Null, None)
            }
            "empty" => (StdLoc::Empty, keep(self.r().seq_id(&variant)?)),
            "whole" => (StdLoc::Whole, keep(self.r().seq_id(&variant)?)),
            "int" => {
                let Interval {
                    from,
                    to,
                    strand,
                    id,
                } = self.interval(&variant)?;
                (StdLoc::Int { from, to, strand }, keep(id))
            }
            "pnt" => {
                const NAMES: [&str; 4] = ["point", "strand", "id", "fuzz"];
                let mut members = self.r().begin_members(&variant, &NAMES)?;
                let (mut point, mut strand, mut id) = (None, None, None);
                while let Some(index) = self.r().next_member(&mut members)? {
                    match index {
                        0 => point = Some(self.r().uint4()?),
                        1 => strand = Some(self.r().enumerated(&NA_STRAND)?.1),
                        2 => id = Some(self.r().seq_id(members.at())?),
                        _ => self.typed_member(named::SEQ_POINT, members.at())?,
                    }
                }
                let (Some(point), Some(id)) = (point, id) else {
                    return Err(self.r().error("a mandatory member is missing"));
                };
                (StdLoc::Pnt { point, strand }, keep(id))
            }
            _ => {
                self.typed_member(named::SEQ_LOC, &variant)?;
                (StdLoc::Other, None)
            }
        };
        Ok((location, id, name))
    }

    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:332-337
    fn messages(&mut self, at: &Frame<'_>) -> Result<()> {
        let mut items = self.r().begin_list(at)?;
        while self.r().next_item(&mut items)? {
            let mut members = self.r().begin_members(items.at(), &["code", "message"])?;
            let mut code = false;
            while let Some(index) = self.r().next_member(&mut members)? {
                if index == 0 {
                    self.r().int4()?;
                    code = true;
                } else {
                    self.string(members.at())?;
                }
            }
            if !code {
                return Err(self.r().error("member code missing"));
            }
        }
        self.r().end_block()
    }

    /// A `Bioseq-set`: its `seq-set`.
    // NCBI reference (598d8ae6): c++/src/objects/seqset/seqset.asn:22-60
    fn bioseq_set(&mut self, at: &Frame<'_>) -> Result<Vec<SeqEntry>> {
        const NAMES: [&str; 9] = [
            "id", "coll", "level", "class", "release", "date", "descr", "seq-set", "annot",
        ];
        let mut members = self.r().begin_members(at, &NAMES)?;
        let mut entries = None;
        while let Some(index) = self.r().next_member(&mut members)? {
            let at = members.at();
            match index {
                0 => {
                    self.r().object_id(at)?;
                }
                1 => {
                    self.r().dbtag(at)?;
                }
                2 => {
                    self.r().int4()?;
                }
                3 => {
                    self.r().enumerated(&BIOSEQ_SET_CLASS)?;
                }
                4 => {
                    self.string(at)?;
                }
                5 => {
                    self.r().date(at)?;
                }
                6 => {
                    let mut ignored = Query::default();
                    self.seq_descr(at, &mut ignored)?;
                }
                7 => {
                    let mut list = Vec::new();
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        let entry = self.r().choice(items.at(), &["seq", "set"])?;
                        list.push(if entry.index() == 0 {
                            SeqEntry::Seq(self.bioseq(&entry)?)
                        } else {
                            self.bioseq_set(&entry)?;
                            SeqEntry::Set
                        });
                    }
                    self.r().end_block()?;
                    entries = Some(list);
                }
                _ => self.typed_member(named::BIOSEQ_SET, at)?,
            }
        }
        entries.ok_or_else(|| self.r().error("member seq-set missing"))
    }

    // NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:36-40
    fn bioseq(&mut self, at: &Frame<'_>) -> Result<Query> {
        let mut query = Query::default();
        let mut members = self
            .r()
            .begin_members(at, &["id", "descr", "inst", "annot"])?;
        let (mut ids, mut inst) = (false, false);
        while let Some(index) = self.r().next_member(&mut members)? {
            let at = members.at();
            match index {
                0 => {
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        let id = self.r().seq_id(items.at())?;
                        query.ids.push(id);
                    }
                    self.r().end_block()?;
                    ids = true;
                }
                1 => self.seq_descr(at, &mut query)?,
                2 => {
                    self.seq_inst(at, &mut query)?;
                    inst = true;
                }
                _ => {
                    query.has_annot = true;
                    self.typed_member(named::BIOSEQ, at)?;
                }
            }
        }
        if !(ids && inst) {
            return Err(self.r().error("a mandatory member is missing"));
        }
        Ok(query)
    }

    /// A `Seq-descr`, its first title and the descriptors LOSAT does not read noted in
    /// `query`. An empty one is NCBI's error.
    // NCBI reference (598d8ae6): c++/src/objects/seq/Seq_descr.cpp:62-68
    // (`OBJECTS_SEQ_DESCR_ALLOW_EMPTY` is false unless the environment sets it, which LOSAT's
    // rpsbproc rejects (`reader_variable`), or a registry file does, which
    // the application start rejects, `ApplicationStart::settings`.)
    fn seq_descr(&mut self, at: &Frame<'_>, query: &mut Query) -> Result<()> {
        let mut count = 0usize;
        let mut items = self.r().begin_list(at)?;
        while self.r().next_item(&mut items)? {
            count += 1;
            let variant = self.r().choice(items.at(), &SEQDESC)?;
            match SEQDESC[variant.index()] {
                "title" => {
                    let title = self.string(&variant)?;
                    if query.title.is_none() {
                        query.title = Some(title);
                    }
                }
                "user" => {
                    if self.user_object(&variant, true)? && query.unread_descriptor.is_none() {
                        query.unread_descriptor =
                            Some("a CFastaReader user object whose DefLine is not a string");
                    }
                }
                _ => self.typed_member(named::SEQDESC, &variant)?,
            }
        }
        self.r().end_block()?;
        if count == 0 {
            return Err(self.r().error("empty Seq-descr is not allowed"));
        }
        Ok(())
    }

    /// A `User-object`: whether it is a `CFastaReader` object with a `DefLine` field whose
    /// value is not a string.
    // NCBI reference (598d8ae6): c++/src/objects/general/general.asn:128-148
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:504-523
    // The fields are read as the module has them; a `str` or `strs` is a UTF8String, which NCBI
    // reads with any byte (`Reader::utf8_string`).
    //
    // `top`: the object is a descriptor's (whose `CFastaReader` DefLine rpsbproc reads). Returns
    // whether it is a `CFastaReader` object with a `DefLine` field whose value is not a string.
    fn user_object(&mut self, at: &Frame<'_>, top: bool) -> Result<bool> {
        let mut members = self.r().begin_members(at, &["class", "type", "data"])?;
        let (mut fasta_reader, mut type_seen, mut data_seen, mut bad_defline) =
            (false, false, false, false);
        while let Some(index) = self.r().next_member(&mut members)? {
            match index {
                0 => {
                    self.string(members.at())?;
                }
                1 => {
                    type_seen = true;
                    let id = self.r().object_id(members.at())?;
                    fasta_reader =
                        top && matches!(id, ObjectId::Str(text) if text == b"CFastaReader");
                }
                _ => {
                    data_seen = true;
                    let mut items = self.r().begin_list(members.at())?;
                    while self.r().next_item(&mut items)? {
                        bad_defline |= self.user_field(items.at(), fasta_reader)?;
                    }
                    self.r().end_block()?;
                }
            }
        }
        if !(type_seen && data_seen) {
            return Err(self.r().error("a mandatory member is missing"));
        }
        Ok(bad_defline)
    }

    /// A `User-field` (`User-object` above); `fasta_reader`: the field is one of a descriptor's
    /// `CFastaReader` object. Returns whether it is a `DefLine` field whose value is not a string.
    fn user_field(&mut self, at: &Frame<'_>, fasta_reader: bool) -> Result<bool> {
        const FIELD_DATA: [&str; 12] = [
            "str", "int", "real", "bool", "os", "object", "strs", "ints", "reals", "oss", "fields",
            "objects",
        ];
        let mut members = self.r().begin_members(at, &["label", "num", "data"])?;
        let (mut defline, mut label, mut data, mut bad_defline) = (false, false, false, false);
        while let Some(index) = self.r().next_member(&mut members)? {
            match index {
                0 => {
                    label = true;
                    let id = self.r().object_id(members.at())?;
                    defline = matches!(id, ObjectId::Str(text) if text == b"DefLine");
                }
                1 => {
                    self.r().int4()?;
                }
                _ => {
                    data = true;
                    let variant = self.r().choice(members.at(), &FIELD_DATA)?;
                    bad_defline = fasta_reader && defline && variant.index() != 0;
                    match FIELD_DATA[variant.index()] {
                        "str" => {
                            self.r().utf8_string()?;
                        }
                        "int" => {
                            self.r().int4()?;
                        }
                        "real" => {
                            self.real()?;
                        }
                        "bool" => {
                            self.r().boolean()?;
                        }
                        "os" => {
                            self.r().octets()?;
                        }
                        "object" => {
                            self.user_object(&variant, false)?;
                        }
                        "strs" => self
                            .r()
                            .each(&variant, |reader, _| reader.utf8_string().map(drop))?,
                        "ints" => self
                            .r()
                            .each(&variant, |reader, _| reader.int4().map(drop))?,
                        "reals" => self
                            .r()
                            .each(&variant, |reader, _| reader.real().map(drop))?,
                        "oss" => self
                            .r()
                            .each(&variant, |reader, _| reader.octets().map(drop))?,
                        "fields" => {
                            let mut items = self.r().begin_list(&variant)?;
                            while self.r().next_item(&mut items)? {
                                self.user_field(items.at(), false)?;
                            }
                            self.r().end_block()?;
                        }
                        _ => {
                            let mut items = self.r().begin_list(&variant)?;
                            while self.r().next_item(&mut items)? {
                                self.user_object(items.at(), false)?;
                            }
                            self.r().end_block()?;
                        }
                    }
                }
            }
        }
        if !(label && data) {
            return Err(self.r().error("a mandatory member is missing"));
        }
        Ok(bad_defline)
    }

    // NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:262-297
    // `length` is `TSeqPos` (unsigned) in NCBI's class.
    fn seq_inst(&mut self, at: &Frame<'_>, query: &mut Query) -> Result<()> {
        const NAMES: [&str; 9] = [
            "repr", "mol", "length", "fuzz", "topology", "strand", "seq-data", "ext", "hist",
        ];
        let mut members = self.r().begin_members(at, &NAMES)?;
        let (mut repr, mut mol) = (false, false);
        while let Some(index) = self.r().next_member(&mut members)? {
            match index {
                0 => {
                    self.r().enumerated(&SEQ_INST_REPR)?;
                    repr = true;
                }
                1 => {
                    query.mol = self.r().enumerated(&SEQ_INST_MOL)?.1;
                    mol = true;
                }
                2 => query.length = Some(self.r().uint4()?),
                4 => {
                    self.r().enumerated(&SEQ_INST_TOPOLOGY)?;
                }
                5 => {
                    self.r().enumerated(&SEQ_INST_STRAND)?;
                }
                6 => query.seq_data = Some(self.seq_data(members.at())?),
                7 => self.seq_ext(members.at(), query)?,
                _ => self.typed_member(named::SEQ_INST, members.at())?,
            }
        }
        if !(repr && mol) {
            return Err(self.r().error("a mandatory member is missing"));
        }
        Ok(())
    }

    /// A `Seq-data`: its variant and its letters (or bytes).
    fn seq_data(&mut self, at: &Frame<'_>) -> Result<(&'static str, Vec<u8>)> {
        let variant = self.r().choice(at, &SEQ_DATA)?;
        let name = SEQ_DATA[variant.index()];
        let data = match name {
            "iupacna" | "iupacaa" | "ncbieaa" => self.string(&variant)?,
            "gap" => {
                self.typed_member(named::SEQ_DATA, &variant)?;
                Vec::new()
            }
            _ => self.r().octets()?,
        };
        Ok((name, data))
    }

    // NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:303-324
    // `Get1LetterSeqData` reads the literals' data of a `delta` (`query::fill_sequence`); the
    // other variants and the `loc` items give no residues and are read as values, and a
    // literal of `gap` data is skipped (`SeqDataTo1LtrSeq` gives it no letters and -1, which
    // `Get1LetterSeqData` does not append). The `length` of a literal is `TSeqPos` (unsigned,
    // seq.def:13-14).
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:151,193-194
    fn seq_ext(&mut self, at: &Frame<'_>, query: &mut Query) -> Result<()> {
        const NAMES: [&str; 4] = ["seg", "ref", "map", "delta"];
        let variant = self.r().choice(at, &NAMES)?;
        if variant.index() != 3 {
            return self.typed_member(named::SEQ_EXT, &variant);
        }
        let mut items = self.r().begin_list(&variant)?;
        while self.r().next_item(&mut items)? {
            let item = self.r().choice(items.at(), &["loc", "literal"])?;
            if item.index() == 0 {
                self.typed_member(named::DELTA_SEQ, &item)?;
                continue;
            }
            let mut members = self
                .r()
                .begin_members(&item, &["length", "fuzz", "seq-data"])?;
            let mut length = false;
            while let Some(index) = self.r().next_member(&mut members)? {
                match index {
                    0 => {
                        self.r().uint4()?;
                        length = true;
                    }
                    1 => self.typed_member(named::SEQ_LITERAL, members.at())?,
                    _ => {
                        let data = self.seq_data(members.at())?;
                        if data.0 != "gap" {
                            query.delta_literals.push(data);
                        }
                    }
                }
            }
            if !length {
                return Err(self.r().error("a mandatory member is missing"));
            }
        }
        self.r().end_block()
    }

    /// A `Seq-align-set` (SET OF Seq-align).
    fn seq_align_set(&mut self, at: &Frame<'_>) -> Result<Vec<Align>> {
        let mut aligns = Vec::new();
        let mut items = self.r().begin_list(at)?;
        while self.r().next_item(&mut items)? {
            aligns.push(self.seq_align(items.at())?);
        }
        self.r().end_block()?;
        Ok(aligns)
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:22-51
    fn seq_align(&mut self, at: &Frame<'_>) -> Result<Align> {
        const NAMES: [&str; 7] = ["type", "dim", "score", "segs", "bounds", "id", "ext"];
        let mut align = Align {
            segs: Segs::Other("dendiag"),
            score: None,
            e_value: None,
            bit_score: None,
            num_ident: None,
            score_of_other_type: None,
        };
        let (mut kind, mut segs) = (false, false);
        let mut members = self.r().begin_members(at, &NAMES)?;
        while let Some(index) = self.r().next_member(&mut members)? {
            let at = members.at();
            match index {
                0 => {
                    self.r().enumerated(&SEQ_ALIGN_TYPE)?;
                    kind = true;
                }
                1 => {
                    self.r().int4()?;
                }
                2 => {
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        self.align_score(items.at(), &mut align)?;
                    }
                    self.r().end_block()?;
                }
                3 => {
                    segs = true;
                    let variant = self.r().choice(at, &SEQ_ALIGN_SEGS)?;
                    align.segs = match variant.index() {
                        1 => Segs::Denseg(self.dense_seg(&variant)?),
                        2 => Segs::Std(Box::new(self.std_segs(&variant)?)),
                        3 => {
                            self.packed_seg(&variant)?;
                            Segs::Other("packed")
                        }
                        index => {
                            self.typed_member(seq_align_segs(), &variant)?;
                            Segs::Other(SEQ_ALIGN_SEGS[index])
                        }
                    };
                }
                5 => {
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        self.r().object_id(items.at())?;
                    }
                    self.r().end_block()?;
                }
                _ => self.typed_member(named::SEQ_ALIGN, at)?,
            }
        }
        if !(kind && segs) {
            return Err(self.r().error("a mandatory member is missing"));
        }
        Ok(align)
    }

    /// One `Score` of a Seq-align, into the fields `ParseAlignScores` reads.
    // NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:283-289
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:214-262
    fn align_score(&mut self, at: &Frame<'_>, align: &mut Align) -> Result<()> {
        let mut members = self.r().begin_members(at, &["id", "value"])?;
        let mut name: Option<&'static str> = None;
        let mut value_seen = false;
        while let Some(index) = self.r().next_member(&mut members)? {
            if index == 0 {
                name = match self.r().object_id(members.at())? {
                    ObjectId::Str(text) => match text.as_slice() {
                        b"score" => Some("score"),
                        b"e_value" => Some("e_value"),
                        b"bit_score" => Some("bit_score"),
                        b"num_ident" => Some("num_ident"),
                        _ => None,
                    },
                    ObjectId::Id(_) => None,
                };
                continue;
            }
            value_seen = true;
            let is_real = self.r().choice(members.at(), &["real", "int"])?.index() == 0;
            let wants_real = matches!(name, Some("e_value" | "bit_score"));
            match (name, is_real) {
                (None, true) => {
                    self.real()?;
                }
                (None, false) => {
                    self.r().int4()?;
                }
                (Some(name), true) if wants_real => {
                    let value = self.real()?;
                    if name == "e_value" {
                        align.e_value = Some(value);
                    } else {
                        align.bit_score = Some(value);
                    }
                }
                (Some(name), false) if !wants_real => {
                    let value = self.r().int4()?;
                    if name == "score" {
                        align.score = Some(value);
                    } else {
                        align.num_ident = Some(value);
                    }
                }
                (Some(name), is_real) => {
                    if is_real {
                        self.real()?;
                    } else {
                        self.r().int4()?;
                    }
                    if align.score_of_other_type.is_none() {
                        align.score_of_other_type = Some(name);
                    }
                }
            }
        }
        if !value_seen {
            return Err(self.r().error("member value missing"));
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:74-81
    // `starts` are `TSignedSeqPos` and `lens` `TSeqPos` (unsigned) in NCBI's class.
    fn dense_seg(&mut self, at: &Frame<'_>) -> Result<DenseSeg> {
        const NAMES: [&str; 7] = [
            "dim", "numseg", "ids", "starts", "lens", "strands", "scores",
        ];
        let mut seg = DenseSeg {
            ids: Vec::new(),
            starts: Vec::new(),
            lens: Vec::new(),
        };
        let mut seen = [false; 7];
        let mut members = self.r().begin_members(at, &NAMES)?;
        while let Some(index) = self.r().next_member(&mut members)? {
            seen[index] = true;
            let at = members.at();
            match index {
                0 | 1 => {
                    self.r().int4()?;
                }
                2 => {
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        let id = self.r().seq_id(items.at())?;
                        seg.ids.push(id);
                    }
                    self.r().end_block()?;
                }
                3 => {
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        seg.starts.push(self.r().int4()?);
                    }
                    self.r().end_block()?;
                }
                4 => {
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        seg.lens.push(self.r().uint4()?);
                    }
                    self.r().end_block()?;
                }
                5 => {
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        self.r().enumerated(&NA_STRAND)?;
                    }
                    self.r().end_block()?;
                }
                _ => {
                    let mut scratch = Align {
                        segs: Segs::Other("dendiag"),
                        score: None,
                        e_value: None,
                        bit_score: None,
                        num_ident: None,
                        score_of_other_type: None,
                    };
                    let mut items = self.r().begin_list(at)?;
                    while self.r().next_item(&mut items)? {
                        self.align_score(items.at(), &mut scratch)?;
                    }
                    self.r().end_block()?;
                }
            }
        }
        if !(seen[1] && seen[2] && seen[3] && seen[4]) {
            return Err(self.r().error("a mandatory member is missing"));
        }
        Ok(seg)
    }
}

#[cfg(test)]
#[path = "archive_tests.rs"]
mod tests;
