//! The reader of rpsbproc's input: ASN.1 text `Blast4-archive` objects, one after another
//! from one stream (rpsblast `-outfmt 11` writes one per query batch), read one object at a time
//! as NCBI's rpsbproc reads them.
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:780-792
//! ```c++
//!         list<TDomQuery> qobjs;
//!         while (m_istream->good())
//!         {
//!             CBlast4_archive objBlastOutput;
//!             try
//!             {
//!                 ObjStreamIn<CBlast4_archive> (*m_istream, objBlastOutput, eSerial_AsnText);
//!             }
//!             catch (...)
//!             {
//!                 break;
//!             }
//! ```
//!
//! An object is read whole before rpsbproc uses it, and any failure of NCBI's typed reader (a
//! member or variant that the module does not have, a member out of order or twice, a missing
//! mandatory member, a `,` before `}`, a value of the wrong form, the end of the text inside the
//! object, anything after the last object that is not an object) ends the input: the objects
//! before it are processed and the run ends normally. The grammar below follows the ASN.1
//! modules for every type that the reader descends into (names, order, OPTIONAL and DEFAULT
//! members, CHOICE variants, ENUMERATED values). The values that rpsbproc does not use are
//! read as values (`Reader::skip_value`) where their type is large, only inside types whose
//! grammar is checked; inside such a skipped value LOSAT accepts what NCBI's typed reader may
//! refuse (the members and variants of: the request bodies other than `queue-search`, the
//! `Seqdesc` variants other than `title` and `user`, `Seq-annot`, `Int-fuzz`, the `Seq-loc`
//! variants other than `int` and `packed-int`, the `Seq-align` members `bounds` and `ext`,
//! the segments other than `denseg`, the large `Blast4-value` variants, the `phi-alignments`,
//! `pssm` and `simple-results` of the results, the `Seq-ext` variants other than `delta` and
//! the `loc` items of a `delta`, and `Seq-hist` and `gap` data).
//!
//! What LOSAT does not read is noted while the object is read and refused once the object is
//! read whole (`ArchiveStream::next_archive`): NCBI may still fail to read the object, which then
//! ends the input without an error.

use std::io::Read;

use crate::report::serial::asn_text_reader::{
    is_first_id_char, AsnTextError, AsnTextErrorKind, Reader, Result, StreamSource,
};
use crate::report::serial::pssm::{SEQ_INST_MOL, SEQ_INST_REPR};
use crate::seqid::model::{ObjectId, ParsedSeqId};

/// One `Blast4-archive`, with what rpsbproc takes from it.
// NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:67-74
// ```text
// Blast4-archive ::= SEQUENCE {
//     request                 Blast4-request,
//     results                 Blast4-get-search-results-reply,
//     messages                SEQUENCE OF Blast4-error OPTIONAL
// }
// ```
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

/// `Seq-align.segs`: a `Dense-seg`, or the name of another variant.
#[derive(Debug)]
pub(crate) enum Segs {
    Denseg(DenseSeg),
    Other(&'static str),
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

/// The objects of one input stream.
pub(crate) struct ArchiveStream<R: Read> {
    reader: Reader<StreamSource<R>>,
    done: bool,
}

/// The `Blast4-request-body` variants.
// NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:76-98
// ```text
// Blast4-request-body ::= CHOICE {
//     finish-params           Blast4-finish-params-request,
//     get-databases           NULL,
//     get-matrices            NULL,
//     get-parameters          NULL,
//     get-paramsets           NULL,
//     get-programs            NULL,
//     get-search-results      Blast4-get-search-results-request,
//     get-sequences           Blast4-get-sequences-request,
//     queue-search            Blast4-queue-search-request,
//     get-request-info        Blast4-get-request-info-request,
//     get-sequence-parts      Blast4-get-seq-parts-request,
//     get-windowmasked-taxids NULL,
//     get-protocol-info	    Blast4-get-protocol-info-request,
//     get-search-info	    Blast4-get-search-info-request,
//     get-databases-ex           Blast4-get-databases-ex-request
// }
// ```
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
// ```text
// Blast4-value ::= CHOICE {
//     big-integer             BigInt,
//     bioseq                  Bioseq,
//     boolean                 BOOLEAN,
//     cutoff                  Blast4-cutoff,
//     integer                 INTEGER,
//     matrix                  PssmWithParameters,
//     real                    REAL,
//     seq-align               Seq-align,
//     seq-id                  Seq-id,
//     seq-loc                 Seq-loc,
//     strand-type             Blast4-strand-type,
//     string                  VisibleString,
//     big-integer-list        SEQUENCE OF BigInt,
//     bioseq-list             SEQUENCE OF Bioseq,
//     boolean-list            SEQUENCE OF BOOLEAN,
//     cutoff-list             SEQUENCE OF Blast4-cutoff,
//     integer-list            SEQUENCE OF INTEGER,
//     matrix-list             SEQUENCE OF PssmWithParameters,
//     real-list               SEQUENCE OF REAL,
//     seq-align-list          SEQUENCE OF Seq-align,
//     seq-id-list             SEQUENCE OF Seq-id,
//     seq-loc-list            SEQUENCE OF Seq-loc,
//     strand-type-list        SEQUENCE OF Blast4-strand-type,
//     string-list             SEQUENCE OF VisibleString,
//     bioseq-set              Bioseq-set,
//     seq-align-set           Seq-align-set,
//     query-mask              Blast4-mask
// }
// ```
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
// ```text
// Blast4-strand-type ::= ENUMERATED {
//     forward-strand          (1),
//     reverse-strand          (2),
//     both-strands            (3)
// }
// ```
const BLAST4_STRAND_TYPE: [(&str, i32); 3] = [
    ("forward-strand", 1),
    ("reverse-strand", 2),
    ("both-strands", 3),
];

// NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:423-431
// ```text
// Blast4-frame-type ::= ENUMERATED {
//     notset                  (0),
//     plus1                   (1),
//     plus2                   (2),
//     plus3                   (3),
//     minus1                  (4),
//     minus2                  (5),
//     minus3                  (6)
// }
// ```
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
// ```text
// Na-strand ::= ENUMERATED {          -- strand of nucleic acid
//     unknown (0) ,
//     plus (1) ,
//     minus (2) ,
//     both (3) ,                -- in forward orientation
//     both-rev (4) ,            -- in reverse orientation
//     other (255) }
// ```
const NA_STRAND: [(&str, i32); 6] = [
    ("unknown", 0),
    ("plus", 1),
    ("minus", 2),
    ("both", 3),
    ("both-rev", 4),
    ("other", 255),
];

// NCBI reference (598d8ae6): c++/src/objects/seqset/seqset.asn:26-51
// ```text
//     class ENUMERATED {
//         not-set (0) ,
//         nuc-prot (1) ,              -- nuc acid and coded proteins
//         ...
//         small-genome-set (23) ,     -- viral segments or mitochondrial minicircles
//         other (255) } DEFAULT not-set ,
// ```
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
// ```text
//     topology ENUMERATED {          -- topology of molecule
//         not-set (0) ,
//         linear (1) ,
//         circular (2) ,
//         tandem (3) ,               -- some part of circle
//         other (255) } DEFAULT linear ,
//     strand ENUMERATED {            -- strandedness in living organism
//         not-set (0) ,
//         ss (1) ,                   -- single strand
//         ds (2) ,                   -- double strand
//         mixed (3) ,
//         other (255) } OPTIONAL ,   -- default ds for DNA, ss for RNA, pept
// ```
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
// ```text
// Seqdesc ::= CHOICE {
//     mol-type GIBB-mol ,          -- type of molecule
//     modif SET OF GIBB-mod ,             -- modifiers
//     method GIBB-method ,         -- sequencing method
//     name VisibleString ,         -- a name for this sequence
//     title VisibleString ,        -- a title for this sequence
//     org Org-ref ,                -- if all from one organism
//     comment VisibleString ,      -- a more extensive comment
//     num Numbering ,              -- a numbering system
//     maploc Dbtag ,               -- map location of this sequence
//     pir PIR-block ,              -- PIR specific info
//     genbank GB-block ,           -- GenBank specific info
//     pub Pubdesc ,                -- a reference to the publication
//     region VisibleString ,       -- overall region (globin locus)
//     user User-object ,           -- user defined object
//     sp SP-block ,                -- SWISSPROT specific info
//     dbxref Dbtag ,               -- xref to other databases
//     embl EMBL-block ,            -- EMBL specific information
//     create-date Date ,           -- date entry first created/released
//     update-date Date ,           -- date of last update
//     prf PRF-block ,              -- PRF specific information
//     pdb PDB-block ,              -- PDB specific information
//     het Heterogen ,              -- cofactor, etc associated but not bound
//     source BioSource ,           -- source of materials, includes Org-ref
//     molinfo MolInfo ,            -- info on the molecule and techniques
//     modelev ModelEvidenceSupport -- model evidence for XM records
// }
// ```
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
// ```text
// Seq-data ::= CHOICE {              -- sequence representations
//     iupacna IUPACna ,              -- IUPAC 1 letter nuc acid code
//     iupacaa IUPACaa ,              -- IUPAC 1 letter amino acid code
//     ncbi2na NCBI2na ,              -- 2 bit nucleic acid code
//     ncbi4na NCBI4na ,              -- 4 bit nucleic acid code
//     ncbi8na NCBI8na ,              -- 8 bit extended nucleic acid code
//     ncbipna NCBIpna ,              -- nucleic acid probabilities
//     ncbi8aa NCBI8aa ,              -- 8 bit extended amino acid codes
//     ncbieaa NCBIeaa ,              -- extended ASCII 1 letter aa codes
//     ncbipaa NCBIpaa ,              -- amino acid probabilities
//     ncbistdaa NCBIstdaa,           -- consecutive codes for std aas
//     gap Seq-gap                    -- gap types
// }
// ```
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
// ```text
// Seq-loc ::= CHOICE {
//     null NULL ,           -- not placed
//     empty Seq-id ,        -- to NULL one Seq-id in a collection
//     whole Seq-id ,        -- whole sequence
//     int Seq-interval ,    -- from to
//     packed-int Packed-seqint ,
//     pnt Seq-point ,
//     packed-pnt Packed-seqpnt ,
//     mix Seq-loc-mix ,
//     equiv Seq-loc-equiv ,  -- equivalent sets of locations
//     bond Seq-bond ,
//     feat Feat-id }         -- indirect, through a Seq-feat
// ```
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
// ```text
//     segs CHOICE {                   -- alignment data
//         dendiag SEQUENCE OF Dense-diag ,
//         denseg              Dense-seg ,
//         std     SEQUENCE OF Std-seg ,
//         packed              Packed-seg ,
//         disc                Seq-align-set,
//         spliced             Spliced-seg,
//         sparse              Sparse-seg
//     } ,
// ```
const SEQ_ALIGN_SEGS: [&str; 7] = [
    "dendiag", "denseg", "std", "packed", "disc", "spliced", "sparse",
];

// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:22-29
// ```text
// Seq-align ::= SEQUENCE {
//     type ENUMERATED {
//         not-set (0) ,
//         global (1) ,
//         diags (2) ,     -- unbroken, but not ordered, diagonals
//         partial (3) ,   -- mapping pieces together
//         disc (4) ,      -- discontinuous alignment
//         other (255) } ,
// ```
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
        ArchiveStream {
            reader: Reader::new(StreamSource::new(input, first), "archive", "rpsbproc")
                .refuse_invisible_bytes()
                .defer_refusals(),
            done: false,
        }
    }

    /// The bytes the stream keeps.
    #[cfg(test)]
    pub(crate) fn buffered(&self) -> usize {
        self.reader.src.buffered()
    }

    /// The next object: `Ok(None)` when NCBI's reader fails at it (the end of the input), the
    /// object when it reads, or the refusal of what LOSAT does not read in it.
    pub(crate) fn next_archive(&mut self) -> std::result::Result<Option<Archive>, AsnTextError> {
        if self.done {
            return Ok(None);
        }
        self.reader.refusal = None;
        match self.archive_object() {
            Ok(archive) => match self.reader.refusal.take() {
                Some(refusal) => {
                    self.done = true;
                    Err(refusal)
                }
                None => Ok(Some(archive)),
            },
            Err(error) => {
                self.done = true;
                match error.kind {
                    AsnTextErrorKind::Format => Ok(None),
                    AsnTextErrorKind::Unsupported => Err(error),
                }
            }
        }
    }

    fn real(&mut self) -> Result<f64> {
        self.reader.real()
    }

    fn string(&mut self) -> Result<Vec<u8>> {
        self.reader.string()
    }

    fn skip_value(&mut self) -> Result<()> {
        self.reader.skip_value()
    }

    /// `ReadFileHeader` and the type check of the root reader, then the object.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:709-721
    // ```c++
    // string CObjectIStreamAsn::ReadFileHeader()
    // {
    //     CTempString id = ReadTypeId(SkipWhiteSpace());
    //     string s(id);
    //     if ( SkipWhiteSpace() == ':' &&
    //          m_Input.PeekCharNoEOF(1) == ':' &&
    //          m_Input.PeekCharNoEOF(2) == '=' ) {
    //         m_Input.SkipChars(3);
    //     }
    //     else
    //         ThrowError(fFormatError, "'::=' expected");
    //     return s;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:892-898
    // ```c++
    //     string name = ReadFileHeader();
    //     const string& tname = typeInfo->GetName();
    //     if ( !name.empty() && !tname.empty() && name != tname ) {
    //         ThrowError(fFormatError,
    //                    "incompatible type "+name+"<>"+typeInfo->GetName());
    //     }
    // ```
    // A type name in brackets (`ReadTypeId`, objistrasn.cpp:292-313) stops before its `]`, so
    // the `::=` check that follows always fails.
    fn archive_object(&mut self) -> Result<Archive> {
        let r = &mut self.reader;
        let c = r.skip_white_space()?;
        if c == b'[' {
            return Err(r.error("'::=' expected"));
        }
        let id = r.scan_end_of_id(is_first_id_char(c));
        if !id.is_empty() && r.bytes(id) != b"Blast4-archive" {
            return Err(r.error("incompatible type"));
        }
        if r.skip_white_space()? != b':' || r.peek_at(1) != b':' || r.peek_at(2) != b'=' {
            return Err(r.error("'::=' expected"));
        }
        r.pos += 3;
        let mut members = r.begin_members()?;
        let (mut request, mut results) = (None, None);
        while let Some(index) = self
            .reader
            .next_member(&mut members, &["request", "results", "messages"])?
        {
            match index {
                0 => request = Some(self.request()?),
                1 => results = Some(self.results()?),
                _ => self.messages()?,
            }
        }
        let (Some((ident, body)), Some(results)) = (request, results) else {
            return Err(self.reader.error("a mandatory member is missing"));
        };
        Ok(Archive {
            ident,
            body,
            results,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:56-63
    // ```text
    // Blast4-request ::= SEQUENCE {
    //     ident                   VisibleString OPTIONAL,
    //     body                    Blast4-request-body
    // }
    // ```
    fn request(&mut self) -> Result<(Option<Vec<u8>>, RequestBody)> {
        let mut members = self.reader.begin_members()?;
        let (mut ident, mut body) = (None, None);
        while let Some(index) = self.reader.next_member(&mut members, &["ident", "body"])? {
            match index {
                0 => ident = Some(self.string()?),
                _ => {
                    let variant = self.reader.choice(&REQUEST_BODY)?;
                    body = Some(if variant == 8 {
                        RequestBody::QueueSearch(self.queue_search()?)
                    } else {
                        self.skip_value()?;
                        RequestBody::Other
                    });
                }
            }
        }
        let body = body.ok_or_else(|| self.reader.error("member body missing"))?;
        Ok((ident, body))
    }

    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:156-166
    // ```text
    // Blast4-queue-search-request ::= SEQUENCE {
    //     program                 VisibleString,
    //     service                 VisibleString,
    //     queries                 Blast4-queries,
    //     subject                 Blast4-subject,
    //     paramset                VisibleString OPTIONAL,
    //     algorithm-options       Blast4-parameters OPTIONAL,
    //     program-options         Blast4-parameters OPTIONAL,
    //     format-options          Blast4-parameters OPTIONAL
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:138-142
    // ```text
    // Blast4-queries ::= CHOICE {
    //     pssm                    PssmWithParameters,
    //     seq-loc-list            SEQUENCE OF Seq-loc,
    //     bioseq-set              Bioseq-set
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:499-503
    // ```text
    // Blast4-subject ::= CHOICE {
    //     database                VisibleString,
    //     sequences               SEQUENCE OF Bioseq,
    //     seq-loc-list            SEQUENCE OF Seq-loc
    // }
    // ```
    fn queue_search(&mut self) -> Result<QueueSearch> {
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
        let mut members = self.reader.begin_members()?;
        let (mut program, mut service, mut queries, mut subject) = (None, None, None, None);
        let mut search = QueueSearch {
            program: Vec::new(),
            queries: Queries::Other,
            subject: Subject::Other("sequences"),
            has_algorithm_options: false,
            matrix_name: None,
            evalue_threshold: None,
        };
        while let Some(index) = self.reader.next_member(&mut members, &NAMES)? {
            match index {
                0 => program = Some(self.string()?),
                1 => service = Some(self.string()?),
                2 => {
                    const QUERIES: [&str; 3] = ["pssm", "seq-loc-list", "bioseq-set"];
                    let variant = self.reader.choice(&QUERIES)?;
                    queries = Some(if variant == 2 {
                        Queries::BioseqSet(self.bioseq_set()?)
                    } else {
                        self.skip_value()?;
                        Queries::Other
                    });
                }
                3 => {
                    const SUBJECT: [&str; 3] = ["database", "sequences", "seq-loc-list"];
                    let variant = self.reader.choice(&SUBJECT)?;
                    subject = Some(if variant == 0 {
                        Subject::Database(self.string()?)
                    } else {
                        self.skip_value()?;
                        Subject::Other(SUBJECT[variant])
                    });
                }
                4 => {
                    self.string()?;
                }
                5 => {
                    search.has_algorithm_options = true;
                    let (matrix, evalue) = self.parameters()?;
                    search.matrix_name = matrix;
                    search.evalue_threshold = evalue;
                }
                _ => {
                    self.parameters()?;
                }
            }
        }
        let (Some(program), Some(_), Some(queries), Some(subject)) =
            (program, service, queries, subject)
        else {
            return Err(self.reader.error("a mandatory member is missing"));
        };
        search.program = program;
        search.queries = queries;
        search.subject = subject;
        Ok(search)
    }

    /// `Blast4-parameters`: the values of the first parameters named `MatrixName` and
    /// `EvalueThreshold`.
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:464-467,505
    // ```text
    // Blast4-parameter ::= SEQUENCE {
    //     name                    VisibleString,
    //     value                   Blast4-value
    // }
    // ...
    // Blast4-parameters ::= SEQUENCE OF Blast4-parameter
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/blast/Blast4_parameters.cpp:186-198
    // ```c++
    // CRef<CBlast4_parameter>
    // CBlast4_parameters::GetParamByName(const string name) const
    // {
    // 	CRef <CBlast4_parameter> blast4ParamFound;
    //
    // 	ITERATE(Tdata,it,this->Get()) {
    // 		if((*it)->GetName() == name) {
    // 			blast4ParamFound = (*it);
    // 			break;
    // 		}
    // 	}
    // 	return blast4ParamFound;
    // }
    // ```
    fn parameters(&mut self) -> Result<(Option<ParamValue>, Option<ParamValue>)> {
        let (mut matrix, mut evalue) = (None, None);
        let mut block_start = self.reader.start_block()?;
        while self.reader.next_element(&mut block_start)? {
            let mut members = self.reader.begin_members()?;
            let (mut name, mut value) = (None, None);
            while let Some(index) = self.reader.next_member(&mut members, &["name", "value"])? {
                match index {
                    0 => name = Some(self.string()?),
                    _ => value = Some(self.blast4_value()?),
                }
            }
            let (Some(name), Some(value)) = (name, value) else {
                return Err(self.reader.error("a mandatory member is missing"));
            };
            if name == b"MatrixName" && matrix.is_none() {
                matrix = Some(value);
            } else if name == b"EvalueThreshold" && evalue.is_none() {
                evalue = Some(value);
            }
        }
        self.reader.end_block()?;
        Ok((matrix, evalue))
    }

    /// A `Blast4-value` (the variants listed at `BLAST4_VALUE`).
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:375-378
    // ```text
    // Blast4-cutoff ::= CHOICE {
    //     e-value                 REAL,
    //     raw-score               INTEGER
    // }
    // ```
    fn blast4_value(&mut self) -> Result<ParamValue> {
        let variant = self.reader.choice(&BLAST4_VALUE)?;
        let name = BLAST4_VALUE[variant];
        match name {
            "big-integer" => {
                self.reader.int8()?;
            }
            "boolean" => {
                self.reader.boolean()?;
            }
            "cutoff" => {
                return Ok(ParamValue::Cutoff(
                    if self.reader.choice(&["e-value", "raw-score"])? == 0 {
                        Cutoff::EValue(self.real()?)
                    } else {
                        self.reader.int4()?;
                        Cutoff::RawScore
                    },
                ));
            }
            "integer" => {
                self.reader.int4()?;
            }
            "real" => {
                self.real()?;
            }
            "strand-type" => {
                self.reader.enumerated(&BLAST4_STRAND_TYPE)?;
            }
            "string" => return Ok(ParamValue::String(self.string()?)),
            _ => self.skip_value()?,
        }
        Ok(ParamValue::Other(name))
    }

    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:287-299
    // ```text
    // Blast4-get-search-results-reply ::= SEQUENCE {
    //     alignments              Seq-align-set OPTIONAL,
    //     phi-alignments          Blast4-phi-alignments OPTIONAL,
    //     masks                   SEQUENCE OF Blast4-mask OPTIONAL,
    //     ka-blocks               SEQUENCE OF Blast4-ka-block OPTIONAL,
    //     search-stats            SEQUENCE OF VisibleString OPTIONAL,
    //     pssm                    PssmWithParameters OPTIONAL,
    //     simple-results          Blast4-simple-results OPTIONAL
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:433-438,454-457
    // ```text
    // Blast4-ka-block ::= SEQUENCE {
    //     lambda                  REAL,
    //     k                       REAL,
    //     h                       REAL,
    //     gapped                  BOOLEAN
    // }
    // ...
    // Blast4-mask ::= SEQUENCE {
    //     locations               SEQUENCE OF Seq-loc,
    //     frame                   Blast4-frame-type
    // }
    // ```
    fn results(&mut self) -> Result<Results> {
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
        let mut members = self.reader.begin_members()?;
        while let Some(index) = self.reader.next_member(&mut members, &NAMES)? {
            match index {
                0 => results.alignments = Some(self.seq_align_set()?),
                2 => {
                    let mut masks = Vec::new();
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        masks.push(self.blast4_mask()?);
                    }
                    self.reader.end_block()?;
                    results.masks = Some(masks);
                }
                3 => {
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        let mut members = self.reader.begin_members()?;
                        let mut seen = 0;
                        while let Some(index) = self
                            .reader
                            .next_member(&mut members, &["lambda", "k", "h", "gapped"])?
                        {
                            seen += 1;
                            if index == 3 {
                                self.reader.boolean()?;
                            } else {
                                self.real()?;
                            }
                        }
                        if seen != 4 {
                            return Err(self.reader.error("a mandatory member is missing"));
                        }
                    }
                    self.reader.end_block()?;
                }
                4 => {
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        self.string()?;
                    }
                    self.reader.end_block()?;
                }
                _ => self.skip_value()?,
            }
        }
        Ok(results)
    }

    fn blast4_mask(&mut self) -> Result<Mask> {
        let mut members = self.reader.begin_members()?;
        let (mut first, mut locations, mut frame) = (None, false, false);
        while let Some(index) = self
            .reader
            .next_member(&mut members, &["locations", "frame"])?
        {
            if index == 0 {
                locations = true;
                let mut block_start = self.reader.start_block()?;
                while self.reader.next_element(&mut block_start)? {
                    let location = self.mask_location()?;
                    if first.is_none() {
                        first = Some(location);
                    }
                }
                self.reader.end_block()?;
            } else {
                frame = true;
                self.reader.enumerated(&BLAST4_FRAME_TYPE)?;
            }
        }
        if !(locations && frame) {
            return Err(self.reader.error("a mandatory member is missing"));
        }
        Ok(Mask { first })
    }

    /// A `Seq-loc` of a mask: the ids of a `packed-int`, the name of another variant.
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:93-101
    // ```text
    // Seq-interval ::= SEQUENCE {
    //     from INTEGER ,
    //     to INTEGER ,
    //     strand Na-strand OPTIONAL ,
    //     id Seq-id ,    -- WARNING: this used to be optional
    //     fuzz-from Int-fuzz OPTIONAL ,
    //     fuzz-to Int-fuzz OPTIONAL }
    //
    // Packed-seqint ::= SEQUENCE OF Seq-interval
    // ```
    // `from` and `to` are `TSeqPos` (unsigned) in NCBI's classes.
    fn mask_location(&mut self) -> Result<MaskLocation> {
        let variant = self.reader.choice(&SEQ_LOC)?;
        match SEQ_LOC[variant] {
            "int" => {
                self.seq_interval()?;
                Ok(MaskLocation::Other("int"))
            }
            "packed-int" => {
                let mut ids = Vec::new();
                let mut block_start = self.reader.start_block()?;
                while self.reader.next_element(&mut block_start)? {
                    ids.push(self.seq_interval()?);
                }
                self.reader.end_block()?;
                Ok(MaskLocation::PackedInt(ids))
            }
            name => {
                self.skip_value()?;
                Ok(MaskLocation::Other(name))
            }
        }
    }

    /// A `Seq-interval`: its id.
    fn seq_interval(&mut self) -> Result<ParsedSeqId> {
        const NAMES: [&str; 6] = ["from", "to", "strand", "id", "fuzz-from", "fuzz-to"];
        let mut members = self.reader.begin_members()?;
        let (mut from, mut to, mut id) = (false, false, None);
        while let Some(index) = self.reader.next_member(&mut members, &NAMES)? {
            match index {
                0 => {
                    self.reader.uint4()?;
                    from = true;
                }
                1 => {
                    self.reader.uint4()?;
                    to = true;
                }
                2 => {
                    self.reader.enumerated(&NA_STRAND)?;
                }
                3 => id = Some(self.reader.seq_id()?),
                _ => self.skip_value()?,
            }
        }
        match (from, to, id) {
            (true, true, Some(id)) => Ok(id),
            _ => Err(self.reader.error("a mandatory member is missing")),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:332-337
    // ```text
    // Blast4-error ::= SEQUENCE {
    //     code                    INTEGER,
    //     message                 VisibleString OPTIONAL
    // }
    // ```
    fn messages(&mut self) -> Result<()> {
        let mut block_start = self.reader.start_block()?;
        while self.reader.next_element(&mut block_start)? {
            let mut members = self.reader.begin_members()?;
            let mut code = false;
            while let Some(index) = self
                .reader
                .next_member(&mut members, &["code", "message"])?
            {
                if index == 0 {
                    self.reader.int4()?;
                    code = true;
                } else {
                    self.string()?;
                }
            }
            if !code {
                return Err(self.reader.error("member code missing"));
            }
        }
        self.reader.end_block()
    }

    /// A `Bioseq-set`: its `seq-set`.
    // NCBI reference (598d8ae6): c++/src/objects/seqset/seqset.asn:22-60
    // ```text
    // Bioseq-set ::= SEQUENCE {      -- just a collection
    //     id Object-id OPTIONAL ,
    //     coll Dbtag OPTIONAL ,          -- to identify a collection
    //     level INTEGER OPTIONAL ,       -- nesting level
    //     class ENUMERATED {
    //         ...
    //         other (255) } DEFAULT not-set ,
    //     release VisibleString OPTIONAL ,
    //     date Date OPTIONAL ,
    //     descr Seq-descr OPTIONAL ,
    //     seq-set SEQUENCE OF Seq-entry ,
    //     annot SET OF Seq-annot OPTIONAL }
    //
    // Seq-entry ::= CHOICE {
    //         seq Bioseq ,
    //         set Bioseq-set }
    // ```
    fn bioseq_set(&mut self) -> Result<Vec<SeqEntry>> {
        const NAMES: [&str; 9] = [
            "id", "coll", "level", "class", "release", "date", "descr", "seq-set", "annot",
        ];
        let mut members = self.reader.begin_members()?;
        let mut entries = None;
        while let Some(index) = self.reader.next_member(&mut members, &NAMES)? {
            match index {
                0 => {
                    self.reader.object_id()?;
                }
                1 => {
                    self.reader.dbtag()?;
                }
                2 => {
                    self.reader.int4()?;
                }
                3 => {
                    self.reader.enumerated(&BIOSEQ_SET_CLASS)?;
                }
                4 => {
                    self.string()?;
                }
                5 => {
                    self.reader.date()?;
                }
                6 => {
                    let mut ignored = Query::default();
                    self.seq_descr(&mut ignored)?;
                }
                7 => {
                    let mut list = Vec::new();
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        list.push(if self.reader.choice(&["seq", "set"])? == 0 {
                            SeqEntry::Seq(self.bioseq()?)
                        } else {
                            self.bioseq_set()?;
                            SeqEntry::Set
                        });
                    }
                    self.reader.end_block()?;
                    entries = Some(list);
                }
                _ => self.skip_value()?,
            }
        }
        entries.ok_or_else(|| self.reader.error("member seq-set missing"))
    }

    // NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:36-40
    // ```text
    // Bioseq ::= SEQUENCE {
    //     id SET OF Seq-id ,            -- equivalent identifiers
    //     descr Seq-descr OPTIONAL , -- descriptors
    //     inst Seq-inst ,            -- the sequence data
    //     annot SET OF Seq-annot OPTIONAL }
    // ```
    fn bioseq(&mut self) -> Result<Query> {
        let mut query = Query::default();
        let mut members = self.reader.begin_members()?;
        let (mut ids, mut inst) = (false, false);
        while let Some(index) = self
            .reader
            .next_member(&mut members, &["id", "descr", "inst", "annot"])?
        {
            match index {
                0 => {
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        let id = self.reader.seq_id()?;
                        query.ids.push(id);
                    }
                    self.reader.end_block()?;
                    ids = true;
                }
                1 => self.seq_descr(&mut query)?,
                2 => {
                    self.seq_inst(&mut query)?;
                    inst = true;
                }
                _ => {
                    query.has_annot = true;
                    self.skip_value()?;
                }
            }
        }
        if !(ids && inst) {
            return Err(self.reader.error("a mandatory member is missing"));
        }
        Ok(query)
    }

    /// A `Seq-descr`, its first title and the descriptors LOSAT does not read noted in
    /// `query`. An empty one is NCBI's error.
    // NCBI reference (598d8ae6): c++/src/objects/seq/Seq_descr.cpp:62-68
    // ```c++
    // void CSeq_descr::PostRead(void) const
    // {
    //     if ( !s_SeqDescrAllowEmpty->Get() && Get().empty() ) {
    //         NCBI_THROW(CSerialException, eInvalidData,
    //                    "empty Seq-descr is not allowed");
    // ```
    // (`OBJECTS_SEQ_DESCR_ALLOW_EMPTY` is false unless the environment sets it, which LOSAT's
    // rpsbproc rejects (`reader_variable`), or a registry file does, which
    // `check_ncbi_application_settings` rejects.)
    fn seq_descr(&mut self, query: &mut Query) -> Result<()> {
        let mut count = 0usize;
        let mut block_start = self.reader.start_block()?;
        while self.reader.next_element(&mut block_start)? {
            count += 1;
            let variant = self.reader.choice(&SEQDESC)?;
            match SEQDESC[variant] {
                "title" => {
                    let title = self.string()?;
                    if query.title.is_none() {
                        query.title = Some(title);
                    }
                }
                "user" => {
                    if self.user_object(true)? && query.unread_descriptor.is_none() {
                        query.unread_descriptor =
                            Some("a CFastaReader user object whose DefLine is not a string");
                    }
                }
                _ => self.skip_value()?,
            }
        }
        self.reader.end_block()?;
        if count == 0 {
            return Err(self.reader.error("empty Seq-descr is not allowed"));
        }
        Ok(())
    }

    /// A `User-object`: whether it is a `CFastaReader` object with a `DefLine` field whose
    /// value is not a string.
    // NCBI reference (598d8ae6): c++/src/objects/general/general.asn:128-148
    // ```text
    // User-object ::= SEQUENCE {
    //     class VisibleString OPTIONAL ,   -- endeavor which designed this object
    //     type Object-id ,                 -- type of object within class
    //     data SEQUENCE OF User-field }    -- the object itself
    //
    // User-field ::= SEQUENCE {
    //     label Object-id ,                -- field label
    //     num INTEGER OPTIONAL ,           -- required for strs, ints, reals, oss
    //     data CHOICE {                    -- field contents
    //         str UTF8String ,
    //         int INTEGER ,
    //         real REAL ,
    //         bool BOOLEAN ,
    //         os OCTET STRING ,
    //         object User-object ,         -- for using other definitions
    //         strs SEQUENCE OF UTF8String ,
    //         ints SEQUENCE OF INTEGER ,
    //         reals SEQUENCE OF REAL ,
    //         oss SEQUENCE OF OCTET STRING ,
    //         fields SEQUENCE OF User-field ,
    //         objects SEQUENCE OF User-object } }
    // ```
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:504-523
    // ```c++
    // 			if ((*iterSeqDesc)->Which() == CSeqdesc::e_User)
    // 			{
    // 				const CUser_object & rUsrDescr = (*iterSeqDesc)->GetUser();
    // 				if (rUsrDescr.CanGetType())
    // 				{
    // 					const CObject_id& typeId = rUsrDescr.GetType();
    // 					if (typeId.Which() == CObject_id::e_Str && NStr::EqualCase(typeId.GetStr(), "CFastaReader"))
    // 					{
    // 						...
    // 							const CObject_id& label = (*iterUsrField)->GetLabel();
    // 							if (label.Which() == CObject_id::e_Str && NStr::EqualCase(label.GetStr(), "DefLine"))
    // 							{
    // 								const CUser_field::TData& fieldValue = (*iterUsrField)->GetData();
    // 								const string& fieldStrVal = fieldValue.GetStr();
    // ```
    // The fields are read as the module has them; a `str` or `strs` is a UTF8String, which NCBI
    // reads with any byte (`Reader::utf8_string`).
    //
    // `top`: the object is a descriptor's (whose `CFastaReader` DefLine rpsbproc reads). Returns
    // whether it is a `CFastaReader` object with a `DefLine` field whose value is not a string.
    fn user_object(&mut self, top: bool) -> Result<bool> {
        let mut members = self.reader.begin_members()?;
        let (mut fasta_reader, mut type_seen, mut data_seen, mut bad_defline) =
            (false, false, false, false);
        while let Some(index) = self
            .reader
            .next_member(&mut members, &["class", "type", "data"])?
        {
            match index {
                0 => {
                    self.string()?;
                }
                1 => {
                    type_seen = true;
                    let id = self.reader.object_id()?;
                    fasta_reader =
                        top && matches!(id, ObjectId::Str(text) if text == b"CFastaReader");
                }
                _ => {
                    data_seen = true;
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        bad_defline |= self.user_field(fasta_reader)?;
                    }
                    self.reader.end_block()?;
                }
            }
        }
        if !(type_seen && data_seen) {
            return Err(self.reader.error("a mandatory member is missing"));
        }
        Ok(bad_defline)
    }

    /// A `User-field` (`User-object` above); `fasta_reader`: the field is one of a descriptor's
    /// `CFastaReader` object. Returns whether it is a `DefLine` field whose value is not a string.
    fn user_field(&mut self, fasta_reader: bool) -> Result<bool> {
        const FIELD_DATA: [&str; 12] = [
            "str", "int", "real", "bool", "os", "object", "strs", "ints", "reals", "oss", "fields",
            "objects",
        ];
        let mut members = self.reader.begin_members()?;
        let (mut defline, mut label, mut data, mut bad_defline) = (false, false, false, false);
        while let Some(index) = self
            .reader
            .next_member(&mut members, &["label", "num", "data"])?
        {
            match index {
                0 => {
                    label = true;
                    defline = matches!(self.reader.object_id()?, ObjectId::Str(text) if text == b"DefLine");
                }
                1 => {
                    self.reader.int4()?;
                }
                _ => {
                    data = true;
                    let variant = self.reader.choice(&FIELD_DATA)?;
                    bad_defline = fasta_reader && defline && variant != 0;
                    match FIELD_DATA[variant] {
                        "str" => {
                            self.reader.utf8_string()?;
                        }
                        "int" => {
                            self.reader.int4()?;
                        }
                        "real" => {
                            self.real()?;
                        }
                        "bool" => {
                            self.reader.boolean()?;
                        }
                        "os" => {
                            self.reader.octets()?;
                        }
                        "object" => {
                            self.user_object(false)?;
                        }
                        "strs" => self.reader.each(|reader| reader.utf8_string().map(drop))?,
                        "ints" => self.reader.each(|reader| reader.int4().map(drop))?,
                        "reals" => self.reader.each(|reader| reader.real().map(drop))?,
                        "oss" => self.reader.each(|reader| reader.octets().map(drop))?,
                        "fields" => {
                            let mut block_start = self.reader.start_block()?;
                            while self.reader.next_element(&mut block_start)? {
                                self.user_field(false)?;
                            }
                            self.reader.end_block()?;
                        }
                        _ => {
                            let mut block_start = self.reader.start_block()?;
                            while self.reader.next_element(&mut block_start)? {
                                self.user_object(false)?;
                            }
                            self.reader.end_block()?;
                        }
                    }
                }
            }
        }
        if !(label && data) {
            return Err(self.reader.error("a mandatory member is missing"));
        }
        Ok(bad_defline)
    }

    // NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:262-297
    // ```text
    // Seq-inst ::= SEQUENCE {            -- the sequence data itself
    //     repr ENUMERATED {              -- representation class
    //         ...
    //     mol ENUMERATED {               -- molecule class in living organism
    //         ...
    //     length INTEGER OPTIONAL ,      -- length of sequence in residues
    //     fuzz Int-fuzz OPTIONAL ,       -- length uncertainty
    //     topology ENUMERATED {          -- topology of molecule
    //         ...
    //     strand ENUMERATED {            -- strandedness in living organism
    //         ...
    //     seq-data Seq-data OPTIONAL ,   -- the sequence
    //     ext Seq-ext OPTIONAL ,         -- extensions for special types
    //     hist Seq-hist OPTIONAL }       -- sequence history
    // ```
    // `length` is `TSeqPos` (unsigned) in NCBI's class.
    fn seq_inst(&mut self, query: &mut Query) -> Result<()> {
        const NAMES: [&str; 9] = [
            "repr", "mol", "length", "fuzz", "topology", "strand", "seq-data", "ext", "hist",
        ];
        let mut members = self.reader.begin_members()?;
        let (mut repr, mut mol) = (false, false);
        while let Some(index) = self.reader.next_member(&mut members, &NAMES)? {
            match index {
                0 => {
                    self.reader.enumerated(&SEQ_INST_REPR)?;
                    repr = true;
                }
                1 => {
                    query.mol = self.reader.enumerated(&SEQ_INST_MOL)?.1;
                    mol = true;
                }
                2 => query.length = Some(self.reader.uint4()?),
                4 => {
                    self.reader.enumerated(&SEQ_INST_TOPOLOGY)?;
                }
                5 => {
                    self.reader.enumerated(&SEQ_INST_STRAND)?;
                }
                6 => query.seq_data = Some(self.seq_data()?),
                7 => self.seq_ext(query)?,
                _ => self.skip_value()?,
            }
        }
        if !(repr && mol) {
            return Err(self.reader.error("a mandatory member is missing"));
        }
        Ok(())
    }

    /// A `Seq-data`: its variant and its letters (or bytes).
    fn seq_data(&mut self) -> Result<(&'static str, Vec<u8>)> {
        let variant = self.reader.choice(&SEQ_DATA)?;
        let name = SEQ_DATA[variant];
        let data = match name {
            "iupacna" | "iupacaa" | "ncbieaa" => self.string()?,
            "gap" => {
                self.skip_value()?;
                Vec::new()
            }
            _ => self.reader.octets()?,
        };
        Ok((name, data))
    }

    // NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:303-324
    // ```text
    // Seq-ext ::= CHOICE {
    //     seg Seg-ext ,        -- segmented sequences
    //     ref Ref-ext ,        -- hot link to another sequence (a view)
    //     map Map-ext ,        -- ordered map of markers
    //     delta Delta-ext }
    // ...
    // Delta-ext ::= SEQUENCE OF Delta-seq
    //
    // Delta-seq ::= CHOICE {
    //     loc Seq-loc ,       -- point to a sequence
    //     literal Seq-literal }   -- a piece of sequence
    //
    // Seq-literal ::= SEQUENCE {
    //     length INTEGER ,         -- must give a length in residues
    //     fuzz Int-fuzz OPTIONAL , -- could be unsure
    //     seq-data Seq-data OPTIONAL } -- may have the data
    // ```
    // `Get1LetterSeqData` reads the literals' data of a `delta` (`query::fill_sequence`); the
    // other variants and the `loc` items give no residues and are read as values, and a
    // literal of `gap` data is skipped (`SeqDataTo1LtrSeq` gives it no letters and -1, which
    // `Get1LetterSeqData` does not append). The `length` of a literal is `TSeqPos` (unsigned,
    // seq.def:13-14).
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:151,193-194
    // ```c++
    // 	dst.clear();
    // ...
    // 	default:
    // 		iDataType = -1;
    // ```
    fn seq_ext(&mut self, query: &mut Query) -> Result<()> {
        const NAMES: [&str; 4] = ["seg", "ref", "map", "delta"];
        if self.reader.choice(&NAMES)? != 3 {
            return self.skip_value();
        }
        let mut block_start = self.reader.start_block()?;
        while self.reader.next_element(&mut block_start)? {
            if self.reader.choice(&["loc", "literal"])? == 0 {
                self.skip_value()?;
                continue;
            }
            let mut members = self.reader.begin_members()?;
            let mut length = false;
            while let Some(index) = self
                .reader
                .next_member(&mut members, &["length", "fuzz", "seq-data"])?
            {
                match index {
                    0 => {
                        self.reader.uint4()?;
                        length = true;
                    }
                    1 => self.skip_value()?,
                    _ => {
                        let data = self.seq_data()?;
                        if data.0 != "gap" {
                            query.delta_literals.push(data);
                        }
                    }
                }
            }
            if !length {
                return Err(self.reader.error("a mandatory member is missing"));
            }
        }
        self.reader.end_block()
    }

    /// A `Seq-align-set` (SET OF Seq-align).
    fn seq_align_set(&mut self) -> Result<Vec<Align>> {
        let mut aligns = Vec::new();
        let mut block_start = self.reader.start_block()?;
        while self.reader.next_element(&mut block_start)? {
            aligns.push(self.seq_align()?);
        }
        self.reader.end_block()?;
        Ok(aligns)
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:22-51
    // ```text
    // Seq-align ::= SEQUENCE {
    //     type ENUMERATED {
    //         ...
    //         other (255) } ,
    //     dim INTEGER OPTIONAL ,     -- dimensionality
    //     score SET OF Score OPTIONAL ,   -- for whole alignment
    //     segs CHOICE {                   -- alignment data
    //         ...
    //     } ,
    //     bounds SET OF Seq-loc OPTIONAL,
    //     id SEQUENCE OF Object-id OPTIONAL,
    //     ext SEQUENCE OF User-object OPTIONAL
    // }
    // ```
    fn seq_align(&mut self) -> Result<Align> {
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
        let mut members = self.reader.begin_members()?;
        while let Some(index) = self.reader.next_member(&mut members, &NAMES)? {
            match index {
                0 => {
                    self.reader.enumerated(&SEQ_ALIGN_TYPE)?;
                    kind = true;
                }
                1 => {
                    self.reader.int4()?;
                }
                2 => {
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        self.align_score(&mut align)?;
                    }
                    self.reader.end_block()?;
                }
                3 => {
                    segs = true;
                    let variant = self.reader.choice(&SEQ_ALIGN_SEGS)?;
                    align.segs = if variant == 1 {
                        Segs::Denseg(self.dense_seg()?)
                    } else {
                        self.skip_value()?;
                        Segs::Other(SEQ_ALIGN_SEGS[variant])
                    };
                }
                5 => {
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        self.reader.object_id()?;
                    }
                    self.reader.end_block()?;
                }
                _ => self.skip_value()?,
            }
        }
        if !(kind && segs) {
            return Err(self.reader.error("a mandatory member is missing"));
        }
        Ok(align)
    }

    /// One `Score` of a Seq-align, into the fields `ParseAlignScores` reads.
    // NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:283-289
    // ```text
    // Score ::= SEQUENCE {
    //     id Object-id OPTIONAL ,
    //     value CHOICE {
    //         real REAL ,
    //         int INTEGER
    //     }
    // }
    // ```
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:214-262
    // ```c++
    // 			if (rScore.IsSetId() && rScore.CanGetId())
    // 			{
    // 				const CScore::TId& rObjectId = rScore.GetId();  //CObject_id
    //
    // 				if (CObject_id::e_Str == rObjectId.Which())
    // 				{
    // 					const CObject_id::TStr& rStr = rObjectId.GetStr();    //string
    // 					if (rStr == "score")
    // 					{
    // 						if (rScore.IsSetValue() && rScore.CanGetValue())
    // 						{
    // 							dst.m_iScore = rScore.GetValue().GetInt();
    // 						}
    // 					}
    // 					else if (rStr == "e_value")
    // 					{
    // 						...
    // 							dst.m_dEValue = rScore.GetValue().GetReal();
    // 					...
    // 					else if (rStr == "bit_score")
    // 					...
    // 							dst.m_dBitScore = rScore.GetValue().GetReal();
    // 					...
    // 					else if (rStr == "num_ident")
    // 					...
    // 							dst.m_iNumIdent = rScore.GetValue().GetInt();
    // ```
    fn align_score(&mut self, align: &mut Align) -> Result<()> {
        let mut members = self.reader.begin_members()?;
        let mut name: Option<&'static str> = None;
        let mut value_seen = false;
        while let Some(index) = self.reader.next_member(&mut members, &["id", "value"])? {
            if index == 0 {
                name = match self.reader.object_id()? {
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
            let is_real = self.reader.choice(&["real", "int"])? == 0;
            let wants_real = matches!(name, Some("e_value" | "bit_score"));
            match (name, is_real) {
                (None, true) => {
                    self.real()?;
                }
                (None, false) => {
                    self.reader.int4()?;
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
                    let value = self.reader.int4()?;
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
                        self.reader.int4()?;
                    }
                    if align.score_of_other_type.is_none() {
                        align.score_of_other_type = Some(name);
                    }
                }
            }
        }
        if !value_seen {
            return Err(self.reader.error("member value missing"));
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:74-81
    // ```text
    // Dense-seg ::= SEQUENCE {          -- for (multiway) global or partial alignments
    //     dim INTEGER DEFAULT 2 ,       -- dimensionality
    //     numseg INTEGER ,              -- number of segments here
    //     ids SEQUENCE OF Seq-id ,      -- sequences in order
    //     starts SEQUENCE OF INTEGER ,  -- start OFFSETS in ids order within segs
    //     lens SEQUENCE OF INTEGER ,    -- lengths in ids order within segs
    //     strands SEQUENCE OF Na-strand OPTIONAL ,
    //     scores SEQUENCE OF Score OPTIONAL }  -- score for each seg
    // ```
    // `starts` are `TSignedSeqPos` and `lens` `TSeqPos` (unsigned) in NCBI's class.
    fn dense_seg(&mut self) -> Result<DenseSeg> {
        const NAMES: [&str; 7] = [
            "dim", "numseg", "ids", "starts", "lens", "strands", "scores",
        ];
        let mut seg = DenseSeg {
            ids: Vec::new(),
            starts: Vec::new(),
            lens: Vec::new(),
        };
        let mut seen = [false; 7];
        let mut members = self.reader.begin_members()?;
        while let Some(index) = self.reader.next_member(&mut members, &NAMES)? {
            seen[index] = true;
            match index {
                0 | 1 => {
                    self.reader.int4()?;
                }
                2 => {
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        let id = self.reader.seq_id()?;
                        seg.ids.push(id);
                    }
                    self.reader.end_block()?;
                }
                3 => {
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        seg.starts.push(self.reader.int4()?);
                    }
                    self.reader.end_block()?;
                }
                4 => {
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        seg.lens.push(self.reader.uint4()?);
                    }
                    self.reader.end_block()?;
                }
                5 => {
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        self.reader.enumerated(&NA_STRAND)?;
                    }
                    self.reader.end_block()?;
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
                    let mut block_start = self.reader.start_block()?;
                    while self.reader.next_element(&mut block_start)? {
                        self.align_score(&mut scratch)?;
                    }
                    self.reader.end_block()?;
                }
            }
        }
        if !(seen[1] && seen[2] && seen[3] && seen[4]) {
            return Err(self.reader.error("a mandatory member is missing"));
        }
        Ok(seg)
    }
}

#[cfg(test)]
#[path = "archive_tests.rs"]
mod tests;
