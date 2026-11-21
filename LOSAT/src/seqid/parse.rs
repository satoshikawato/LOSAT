//! The text parsers of `CSeq_id` (`c++/src/objects/seqloc/Seq_id.cpp`, `Textseq_id.cpp`,
//! `c++/src/objects/general/Object_id.cpp`) that `makeblastdb -parse_seqids` runs on a
//! definition line: `ParseIDs` (FASTA ids such as `gi|123|ref|NP_000549.1|`), the
//! constructor from a bare id (`CSeq_id(text, flags)`, with the accession guide) and the
//! typed `Set`. Each gives a `ParsedSeqId` or the `CSeqIdException` NCBI throws
//! (`SeqIdError`, whose `what` is the text NCBI prints when it posts one).

use super::accguide::{
    acc_type, identify_accession, parse_flags::*, SeqIdWarning, SEQID_DDBJ, SEQID_EMBL,
    SEQID_GENBANK, SEQID_GENERAL, SEQID_GI, SEQID_GIBBMT, SEQID_GIBBSQ, SEQID_GIIM, SEQID_GPIPE,
    SEQID_LOCAL, SEQID_NAMED_ANNOT_TRACK, SEQID_NOT_SET, SEQID_OTHER, SEQID_PATENT, SEQID_PDB,
    SEQID_PIR, SEQID_PRF, SEQID_SWISSPROT, SEQID_TPD, SEQID_TPE, SEQID_TPG,
};
use super::model::{
    Dbtag, GiimportId, IdPat, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId, PdbSeqId,
    TextseqId, SEQ_ID_ALTERNATIVES,
};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str::{is_blank, string_to_int8, string_to_non_negative_int, truncate_spaces};

/// A `CSeqIdException` (`eFormat`): where NCBI throws it and its message.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqIdError {
    /// The source file in NCBI's build tree (`__FILE__`: `src/objects/seq/../seqloc/Seq_id.cpp`).
    pub file: &'static str,
    /// The line of the `NCBI_THROW`.
    pub line: u32,
    /// The function, as `NCBI_CURRENT_FUNCTION` names it.
    pub function: &'static str,
    /// The exception's message, as bytes: an id's text keeps the bytes of NCBI's `std::string`
    /// (a seqid list's word that is not UTF-8: probe R3 of line LIB-B).
    pub message: Vec<u8>,
}

// makeblastdb 2.17.0 (conda build) printed, for `>xyz|ABC|DEF`:
// `T0 "/opt/conda/conda-bld/blast_1754905306566/work/c++/src/objects/seq/../seqloc/Seq_id.cpp",
// line 2663: Error: (CSeqIdException::eFormat) ncbi::objects::CSeq_id::ParseIDs() - Unsupported
// ID type xyz`. LOSAT names the build tree of the official Linux x64 release
// (`ncbi_exception::NCBI_BUILD_ROOT`); the path is not compared (makedb-v5 DECISIONS 11).
impl SeqIdError {
    fn new(
        file: &'static str,
        line: u32,
        function: &'static str,
        message: impl Into<Vec<u8>>,
    ) -> Self {
        SeqIdError {
            file,
            line,
            function,
            message: message.into(),
        }
    }

    /// The exception NCBI throws (`NCBI_THROW(CSeqIdException, eFormat, message)`), with the
    /// message's bytes.
    pub fn exception(&self) -> NcbiException {
        let site = ThrowSite {
            file: self.file,
            line: self.line,
            kind: "CSeqIdException::eFormat",
            function: Some(self.function),
        };
        NcbiException::new(site, self.message.clone())
    }

    /// `CException::what()`: the report of the exception, ending with a newline.
    pub fn what(&self) -> Vec<u8> {
        self.exception().what()
    }

    /// `CException::what()` in NCBI's words (before `crate::wording::losat_bytes`), with the
    /// message's bytes as NCBI's `std::string` holds them: windowmasker posts it for an id of a
    /// list, which may be any bytes.
    pub(crate) fn ncbi_what_bytes(&self) -> Vec<u8> {
        self.exception().ncbi_what()
    }
}

const SEQ_ID_CPP: &str = "src/objects/seq/../seqloc/Seq_id.cpp";
const TEXTSEQ_ID_CPP: &str = "src/objects/seq/../seqloc/Textseq_id.cpp";

/// `CObject_id::SetStrOrId(str)`.
///
/// NCBI reference (598d8ae6): c++/src/objects/general/Object_id.cpp:187-198
pub fn str_or_id(s: &[u8]) -> ObjectId {
    if matches!(s.first(), Some(b'1'..=b'9')) {
        let id = string_to_non_negative_int(s);
        if id > 0 {
            return ObjectId::Id(id);
        }
    }
    ObjectId::Str(s.to_vec())
}

/// `CSeq_id::WhichInverseSeqId(code)`: the type of a FASTA tag (`E_Choice`, 0 when it is
/// not one), compared without case.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:526-568,598-607
/// (The commented-out aliases of the array are left out.)
pub fn which_inverse_seq_id(code: &[u8]) -> u32 {
    const TAGS: [(&[u8], u32); 23] = [
        (b"bbm", SEQID_GIBBMT),
        (b"bbs", SEQID_GIBBSQ),
        (b"dbj", SEQID_DDBJ),
        (b"emb", SEQID_EMBL),
        (b"gb", SEQID_GENBANK),
        (b"gi", SEQID_GI),
        (b"gibbsq", SEQID_GIBBSQ),
        (b"gim", SEQID_GIIM),
        (b"gnl", SEQID_GENERAL),
        (b"gpp", SEQID_GPIPE),
        (b"lcl", SEQID_LOCAL),
        (b"nat", SEQID_NAMED_ANNOT_TRACK),
        (b"pat", SEQID_PATENT),
        (b"pdb", SEQID_PDB),
        (b"pgp", SEQID_PATENT),
        (b"pir", SEQID_PIR),
        (b"prf", SEQID_PRF),
        (b"ref", SEQID_OTHER),
        (b"sp", SEQID_SWISSPROT),
        (b"tpd", SEQID_TPD),
        (b"tpe", SEQID_TPE),
        (b"tpg", SEQID_TPG),
        (b"tr", SEQID_SWISSPROT),
    ];
    TAGS.iter()
        .find(|(tag, _)| tag.eq_ignore_ascii_case(code))
        .map_or(SEQID_NOT_SET, |&(_, choice)| choice)
}

/// `s_TextId[type]`: the FASTA tag of a type.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:572-597
fn text_id(choice: u32) -> &'static str {
    const TEXT_ID: [&str; 21] = [
        "???", "lcl", "bbs", "bbm", "gim", "gb", "emb", "pir", "sp", "pat", "ref", "gnl", "gi",
        "dbj", "prf", "pdb", "tpg", "tpe", "tpd", "gpp", "nat",
    ];
    TEXT_ID.get(choice as usize).copied().unwrap_or("")
}

/// `CSeq_id::SelectionName(type)`: the ASN.1 name of the alternative.
fn selection_name(choice: u32) -> &'static str {
    match choice {
        0 => "not set",
        _ => SEQ_ID_ALTERNATIVES
            .get(choice as usize - 1)
            .copied()
            .unwrap_or(""),
    }
}

/// `CSeq_id::ETypeVariant`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum TypeVariant {
    Plain,
    Tr,
    Pgp,
}

/// `CSeq_id::x_IdentifyTypeVariant(type, str)`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2438-2449
fn identify_type_variant(choice: u32, s: &[u8]) -> TypeVariant {
    if choice == SEQID_SWISSPROT && s.eq_ignore_ascii_case(b"tr") {
        TypeVariant::Tr
    } else if choice == SEQID_PATENT && s.eq_ignore_ascii_case(b"pgp") {
        TypeVariant::Pgp
    } else {
        TypeVariant::Plain
    }
}

/// `CTextseq_id::Set(acc, name, version, release, allow_dot_version)`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Textseq_id.cpp:58-139
fn textseq_set(
    acc_in: &[u8],
    name_in: &[u8],
    version: i32,
    release_in: &[u8],
    allow_dot_version: bool,
) -> Result<TextseqId, SeqIdError> {
    const FUNCTION: &str = "ncbi::objects::CTextseq_id::Set";
    if version < 0 {
        return Err(SeqIdError::new(
            TEXTSEQ_ID_CPP,
            67,
            FUNCTION,
            [
                format!("Unexpected negative version {version} for accession ").as_bytes(),
                acc_in,
            ]
            .concat(),
        ));
    }
    let acc = truncate_spaces(acc_in);
    let name = truncate_spaces(name_in);
    let release = truncate_spaces(release_in);
    let mut id = TextseqId::default();
    if !acc.is_empty() {
        let idx = if allow_dot_version {
            acc.iter().rposition(|&c| c == b'.')
        } else {
            None
        };
        match idx {
            None => {
                id.accession = Some(acc.to_vec());
                id.version = (version > 0).then_some(version);
            }
            Some(idx) => {
                let ver = string_to_non_negative_int(&acc[idx + 1..]);
                if ver <= 0 {
                    return Err(SeqIdError::new(
                        TEXTSEQ_ID_CPP,
                        108,
                        FUNCTION,
                        [
                            b"Version embedded in accession ".as_slice(),
                            acc,
                            b" is not a positive integer",
                        ]
                        .concat(),
                    ));
                } else if version > 0 && ver != version {
                    return Err(SeqIdError::new(
                        TEXTSEQ_ID_CPP,
                        112,
                        FUNCTION,
                        [
                            format!("Incompatible version {version} supplied for accession ")
                                .as_bytes(),
                            acc,
                        ]
                        .concat(),
                    ));
                }
                id.accession = Some(acc[..idx].to_vec());
                id.version = Some(ver);
            }
        }
    }
    if !name.is_empty() {
        id.name = Some(name.to_vec());
    }
    if acc.is_empty() && name.is_empty() {
        return Err(SeqIdError::new(
            TEXTSEQ_ID_CPP,
            129,
            FUNCTION,
            [
                format!(
                    "Accession and name missing for Textseq-id (but got version {version}, release "
                )
                .as_bytes(),
                release,
                b")",
            ]
            .concat(),
        ));
    } else if !release.is_empty() {
        id.release = Some(release.to_vec());
    }
    Ok(id)
}

/// `CSeq_id::Set(the_type, the_id)` of a numeric id.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2850-2884
/// `the_id` is a non-negative `int` here (`StringToNonNegativeInt`).
fn set_numeric(choice: u32, the_id: i32) -> Result<ParsedSeqId, SeqIdError> {
    const FUNCTION: &str = "ncbi::objects::CSeq_id::Set";
    if the_id < 0 || (choice == SEQID_LOCAL && the_id == 0) {
        return Err(SeqIdError::new(
            SEQ_ID_CPP,
            2854,
            FUNCTION,
            format!("Non-positive numeric ID {the_id}"),
        ));
    }
    Ok(match choice {
        SEQID_LOCAL => ParsedSeqId::Local(ObjectId::Id(the_id)),
        SEQID_GIBBSQ => ParsedSeqId::Gibbsq(the_id),
        SEQID_GIBBMT => ParsedSeqId::Gibbmt(the_id),
        SEQID_GIIM => ParsedSeqId::Giim(GiimportId {
            id: the_id,
            db: None,
            release: None,
        }),
        SEQID_GI => ParsedSeqId::Gi(i64::from(the_id)),
        _ => {
            return Err(SeqIdError::new(
                SEQ_ID_CPP,
                2880,
                FUNCTION,
                format!("Invalid numeric ID type{}", selection_name(choice)),
            ))
        }
    })
}

/// `CSeq_id::Set(the_type, acc, name, version, release)`; `None` for `e_not_set`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2898-3036
/// (`TGi` is an `Int8`: `NCBI_INT8_GI` is defined in 2.17.0. The Info message is not shown
/// by makeblastdb.)
pub fn set_typed(
    choice: u32,
    acc_in: &[u8],
    name_in: &[u8],
    version: i32,
    release_in: &[u8],
) -> Result<Option<ParsedSeqId>, SeqIdError> {
    let acc = truncate_spaces(acc_in);
    let textseq = |allow_dot| textseq_set(acc, name_in, version, release_in, allow_dot);
    Ok(Some(match choice {
        SEQID_NOT_SET => return Ok(None),
        SEQID_LOCAL => ParsedSeqId::Local(str_or_id(acc)),
        SEQID_GIBBSQ | SEQID_GIBBMT | SEQID_GIIM | SEQID_GI => {
            // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2925-2935
            // `StringToNumeric<TGi>` is `NStr::StringToInt8` with no flags; `acc` is a
            // `CTempString` (the accession with its spaces truncated): a NUL byte in it ends the
            // digits and passes the end check.
            if choice == SEQID_GI {
                if let Ok(gi) = string_to_int8(acc, 0) {
                    if gi > 0 {
                        return Ok(Some(ParsedSeqId::Gi(gi)));
                    }
                }
            }
            let the_id = string_to_non_negative_int(acc);
            if the_id >= 0 {
                set_numeric(choice, the_id)?
            } else {
                return Err(SeqIdError::new(
                    SEQ_ID_CPP,
                    2941,
                    "ncbi::objects::CSeq_id::Set",
                    [
                        format!(
                            "Negative, excessively large, or non-numeric {} ID ",
                            selection_name(choice)
                        )
                        .as_bytes(),
                        acc,
                    ]
                    .concat(),
                ));
            }
        }
        SEQID_GENBANK => ParsedSeqId::Genbank(textseq(true)?),
        SEQID_EMBL => ParsedSeqId::Embl(textseq(true)?),
        SEQID_PIR => ParsedSeqId::Pir(textseq(false)?),
        SEQID_SWISSPROT => ParsedSeqId::Swissprot(textseq(true)?),
        SEQID_OTHER => ParsedSeqId::Other(textseq(true)?),
        SEQID_DDBJ => ParsedSeqId::Ddbj(textseq(true)?),
        SEQID_PRF => ParsedSeqId::Prf(textseq(false)?),
        SEQID_TPG => ParsedSeqId::Tpg(textseq(true)?),
        SEQID_TPE => ParsedSeqId::Tpe(textseq(true)?),
        SEQID_TPD => ParsedSeqId::Tpd(textseq(true)?),
        SEQID_GPIPE => ParsedSeqId::Gpipe(textseq(true)?),
        SEQID_NAMED_ANNOT_TRACK => ParsedSeqId::NamedAnnotTrack(textseq(true)?),
        SEQID_PATENT => {
            let name = truncate_spaces(name_in).to_vec();
            let release = truncate_spaces(release_in);
            let id = if release.eq_ignore_ascii_case(b"pgp") {
                PatentNumber::AppNumber(name)
            } else {
                PatentNumber::Number(name)
            };
            ParsedSeqId::Patent(PatentSeqId {
                seqid: version,
                cit: IdPat {
                    country: acc.to_vec(),
                    id,
                    doc_type: None,
                },
            })
        }
        SEQID_GENERAL => ParsedSeqId::General(Dbtag {
            db: acc.to_vec(),
            tag: str_or_id(truncate_spaces(name_in)),
        }),
        SEQID_PDB => {
            let mut name = truncate_spaces(name_in);
            let chain = if name_in.is_empty() {
                None
            } else if name.is_empty() {
                name = b" ";
                Some(i32::from(b' '))
            } else if name.len() == 1 {
                Some(i32::from(name[0]))
            } else {
                None
            };
            ParsedSeqId::Pdb(PdbSeqId {
                mol: acc.to_vec(),
                chain,
                rel: None,
                chain_id: (!name.is_empty()).then(|| name.to_vec()),
            })
        }
        _ => {
            return Err(SeqIdError::new(
                SEQ_ID_CPP,
                3026,
                "ncbi::objects::CSeq_id::Set",
                format!("Unsupported Seq-id type {}", selection_name(choice)),
            ))
        }
    }))
}

/// `s_SplitVersion(acc_in, acc, ver)`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:151-166
fn split_version(acc_in: &[u8]) -> Result<(&[u8], i32), SeqIdError> {
    let Some(dot) = acc_in.iter().position(|&c| c == b'.') else {
        return Ok((acc_in, 0));
    };
    let verstr = &acc_in[dot + 1..];
    if verstr.is_empty() {
        return Ok((&acc_in[..dot], 0));
    }
    let ver = string_to_non_negative_int(verstr);
    if ver <= 0 {
        return Err(SeqIdError::new(
            SEQ_ID_CPP,
            161,
            "ncbi::objects::s_SplitVersion",
            [
                b"Version embedded in accession ".as_slice(),
                acc_in,
                b" is not a positive integer",
            ]
            .concat(),
        ));
    }
    Ok((&acc_in[..dot], ver))
}

/// `NStr::Split(s, "|", pieces)` without flags: every `|` ends a piece, empty ones kept.
fn split_bars(s: &[u8]) -> std::collections::VecDeque<&[u8]> {
    if s.is_empty() {
        return Default::default();
    }
    s.split(|&c| c == b'|').collect()
}

/// `CSeq_id::CheckLocalID(s)` as `IsValidLocalID` reads it: not blank, at most
/// `kMaxLocalIDLength` (50) bytes, no byte that `CSeq_id_find_pred` finds (a control
/// byte, a byte from 0x80 up, or one of ` >[]|"`; [`has_invalid_local_char`]).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2576-2601
pub fn is_valid_local_id(s: &[u8]) -> bool {
    !is_blank(s) && s.len() <= MAX_LOCAL_ID_LENGTH && !has_invalid_local_char(s)
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:77-86
// `char` is signed on the reference platform (x86-64), so every byte from 0x80 up is
// below 32 there; with an unsigned `char` it is above 127. Either way a control byte
// (below 0x20) and a byte of 0x80 or more are illegal; 0x7F is not. `strchr` also finds
// the terminating NUL, which the `ch<32` test has already taken.
/// `CheckLocalID(s) & fInvalidChar`: a byte below 0x20, from 0x80 up, or one of ` >[]|"`.
pub fn has_invalid_local_char(s: &[u8]) -> bool {
    s.iter()
        .any(|&c| !(0x20..0x80).contains(&c) || b" >[]|\"".contains(&c))
}

/// `CSeq_id::kMaxLocalIDLength`, `kMaxGeneralTagLength`, `kMaxAccessionLength`.
///
/// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id.hpp:843-846
pub const MAX_LOCAL_ID_LENGTH: usize = 50;
pub const MAX_GENERAL_TAG_LENGTH: usize = 50;
pub const MAX_ACCESSION_LENGTH: usize = 30;

/// `kSupportedRawDbtags`: the databases of a bare `db:tag` that is a general id.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:96-124
/// The set compares with case and the key is upper-cased first, so `dbGSS` and `dbSTS`
/// never match.
const RAW_DBTAGS: [&[u8]; 25] = [
    b"ATGC",
    b"BCMHGSC",
    b"BERKELEY",
    b"CELERA",
    b"GSDB",
    b"HOOD",
    b"LANLCHGS",
    b"LRG",
    b"MIPS",
    b"NCBI_EXT_ACC",
    b"NCBI_GENOMES",
    b"NCBI_MITO",
    b"PGEC",
    b"PID",
    b"SGD",
    b"SHGC",
    b"SRA",
    b"TIGR",
    b"UOKNOR",
    b"UWGC",
    b"WASHU",
    b"WIBR",
    b"WUGSC",
    b"dbGSS",
    b"dbSTS",
];

/// `s_CheckForFastaTag(s)`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:632-642
fn check_for_fasta_tag(s: &[u8]) -> u32 {
    if s.len() > 3 && s[2] == b'|' {
        which_inverse_seq_id(&s[..2])
    } else if s.len() > 4 && s[3] == b'|' {
        which_inverse_seq_id(&s[..3])
    } else {
        SEQID_NOT_SET
    }
}

/// `s_HasFastaTag(s)` (Seq_id.cpp:610-620): `|` at the third or fourth byte, with more
/// after it.
fn has_fasta_tag(s: &[u8]) -> bool {
    (s.len() > 3 && s[2] == b'|') || (s.len() > 4 && s[3] == b'|')
}

/// `CSeq_id(the_id, flags)` (`CSeq_id::Set(the_id, flags)`). `warn` gets the diagnostics
/// NCBI posts.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2457-2557
/// (`Set(e_Gi, the_id)` is the typed `Set` with the id as the accession. The tagged
/// branch is reached only without `fParse_NoFASTA`, which makeblastdb's two parses never
/// use on a tagged id.)
pub fn seq_id_from_text(
    the_id_in: &[u8],
    flags: u32,
    warn: &mut dyn FnMut(SeqIdWarning),
) -> Result<ParsedSeqId, SeqIdError> {
    const FUNCTION: &str = "ncbi::objects::CSeq_id::Set";
    let the_id = truncate_spaces(the_id_in);
    let mut choice = SEQID_NOT_SET;
    if flags & NO_FASTA == 0 {
        choice = check_for_fasta_tag(the_id);
    }
    let set = |choice, acc: &[u8], name: &[u8], version| {
        set_typed(choice, acc, name, version, b"").map(|id| id.expect("a set type gives an id"))
    };
    if choice == SEQID_NOT_SET {
        if the_id.is_empty() {
            return Err(SeqIdError::new(
                SEQ_ID_CPP,
                2468,
                FUNCTION,
                "Empty bare accession supplied".to_string(),
            ));
        }
        if flags & ANY_RAW != 0 {
            choice = acc_type(identify_accession(the_id, flags | FALLBACK_OK, warn));
        }
        match choice {
            SEQID_GI => set(SEQID_GI, the_id, b"", 0),
            SEQID_NOT_SET => {
                if let Some(colon) = the_id.iter().position(|&c| c == b':') {
                    let db = the_id[..colon].to_ascii_uppercase();
                    if RAW_DBTAGS.contains(&db.as_slice()) {
                        return set(SEQID_GENERAL, &the_id[..colon], &the_id[colon + 1..], 0);
                    }
                }
                if flags & VALID_LOCAL != 0
                    && (flags & ANY_LOCAL == ANY_LOCAL || is_valid_local_id(the_id))
                {
                    set(SEQID_LOCAL, the_id, b"", 0)
                } else {
                    Err(SeqIdError::new(
                        SEQ_ID_CPP,
                        2501,
                        FUNCTION,
                        [b"Malformatted ID ".as_slice(), the_id].concat(),
                    ))
                }
            }
            SEQID_PRF => set(SEQID_PRF, b"", the_id, 0),
            SEQID_PDB => {
                let mol = &the_id[..the_id.len().min(4)];
                let chain: &[u8] = if the_id.len() > 5 {
                    &the_id[5..]
                } else if the_id.len() == 5 && the_id[4] != b'|' {
                    &the_id[4..5]
                } else {
                    b""
                };
                set(SEQID_PDB, mol, chain, 0)
            }
            _ => {
                let (acc, ver) = split_version(the_id)?;
                set(choice, acc, b"", ver)
            }
        }
    } else {
        let mut pieces = split_bars(the_id);
        let tv = identify_type_variant(choice, pieces.pop_front().unwrap_or_default());
        let (id, _) = x_init(&mut pieces, choice, tv)?;
        // One warning for each non-empty piece left (`ITERATE` without `break`); no caller
        // in makeblastdb passes `fParse_PartialOK` here (ParseIDs adds `fParse_NoFASTA`).
        for piece in &pieces {
            if piece.is_empty() {
                continue;
            }
            if flags & PARTIAL_OK != 0 {
                warn(SeqIdWarning {
                    subcode: 10,
                    message: [
                        b"Ignoring extra parts (synonyms?) in FASTA-style ID ".as_slice(),
                        the_id,
                    ]
                    .concat(),
                });
            } else {
                return Err(SeqIdError::new(
                    SEQ_ID_CPP,
                    2542,
                    FUNCTION,
                    [
                        b"FASTA-style ID ".as_slice(),
                        the_id,
                        b" has too many parts.",
                    ]
                    .concat(),
                ));
            }
        }
        Ok(id.expect("a FASTA tag gives a type"))
    }
}

/// `CSeq_id::x_Init(fasta_pieces, type, tv)`: the id of the type from the next pieces,
/// and the type that the piece left first announces (`e_not_set` for none).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2698-2840
/// (`strchr` also finds the NUL of a `fields[0]` whose fifth byte is NUL.)
fn x_init(
    pieces: &mut std::collections::VecDeque<&[u8]>,
    choice: u32,
    tv: TypeVariant,
) -> Result<(Option<ParsedSeqId>, u32), SeqIdError> {
    const FUNCTION: &str = "ncbi::objects::CSeq_id::x_Init";
    let (min_fields, max_fields) = match choice {
        SEQID_LOCAL | SEQID_GIBBSQ | SEQID_GIBBMT | SEQID_GIIM | SEQID_GI => (1, 1),
        SEQID_PATENT => (3, 3),
        SEQID_GENERAL => (2, 2),
        _ => (1, 2),
    };
    let mut fields: [Vec<u8>; 3] = Default::default();
    let mut next_type = SEQID_NOT_SET;
    for i in 0..max_fields {
        let Some(&front) = pieces.front() else {
            if i >= min_fields {
                break;
            }
            return Err(SeqIdError::new(
                SEQ_ID_CPP,
                2738,
                FUNCTION,
                format!("Not enough fields for ID of type {}", text_id(choice)),
            ));
        };
        if i >= min_fields && pieces.len() > 1 && (front.len() == 2 || front.len() == 3) {
            next_type = which_inverse_seq_id(front);
            if next_type != SEQID_NOT_SET {
                let second = pieces[1];
                let next_type_2 = if second.len() == 2 || second.len() == 3 {
                    which_inverse_seq_id(second)
                } else {
                    SEQID_NOT_SET
                };
                if next_type_2 != SEQID_NOT_SET {
                    next_type = next_type_2;
                } else if second.iter().all(|c| b" \t\n".contains(c)) && pieces.len() == 2 {
                    next_type = SEQID_NOT_SET;
                } else {
                    break;
                }
            }
        }
        fields[i] = front.to_vec();
        pieces.pop_front();
    }
    if choice == SEQID_GENERAL && fields[0].eq_ignore_ascii_case(b"dbSNP") && !pieces.is_empty() {
        let mut snp_name = std::mem::take(&mut fields[1]);
        for piece in pieces.drain(..) {
            snp_name.push(b'|');
            snp_name.extend_from_slice(piece);
        }
        fields[1] = snp_name;
    }
    while pieces.front().is_some_and(|piece| piece.is_empty()) {
        pieces.pop_front();
    }
    let mut ver = 0;
    match choice {
        SEQID_SWISSPROT => {
            fields[2] = if tv == TypeVariant::Tr {
                b"unreviewed".to_vec()
            } else {
                b"reviewed".to_vec()
            };
        }
        SEQID_PATENT => {
            ver = string_to_non_negative_int(&fields[2]);
            if ver < 0 {
                return Err(SeqIdError::new(
                    SEQ_ID_CPP,
                    2799,
                    FUNCTION,
                    [
                        b"Bad sequence number ".as_slice(),
                        &fields[2],
                        b" for ",
                        &fields[0],
                        b" patent ",
                        &fields[1],
                    ]
                    .concat(),
                ));
            }
            fields[2] = b"pat".to_vec();
        }
        SEQID_PDB => {
            let mol = &fields[0];
            if mol.len() < 4
                || (mol.len() > 5
                    && (!fields[1].is_empty() || !matches!(mol[4], b'|' | b'-' | b'_' | 0)))
            {
                return Err(SeqIdError::new(
                    SEQ_ID_CPP,
                    2813,
                    FUNCTION,
                    [b"Malformatted PDB ID ".as_slice(), mol].concat(),
                ));
            }
            if mol.len() > 4 && fields[1].is_empty() {
                fields[1] = if mol.len() > 5 {
                    mol[5..].to_vec()
                } else {
                    mol[4..].to_vec()
                };
                fields[0].truncate(4);
            }
        }
        _ => {}
    }
    let id = set_typed(choice, &fields[0], &fields[1], ver, &fields[2])?;
    Ok((id, next_type))
}

/// `CSeq_id::ParseIDs(ids, s, flags)`: the ids of a FASTA id string; an id that stays
/// unset (`None`) after a skipped piece is kept, as NCBI keeps it. `warn` gets the
/// exceptions that `fParse_PartialOK` turns into warnings (subcode 7, the `what` text).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2613-2696
/// A throw of `x_Init` leaves `type` as it was (the assignment does not happen), so the
/// next pieces are read as that type again.
pub fn parse_ids(
    s: &[u8],
    flags: u32,
    warn: &mut dyn FnMut(SeqIdWarning),
) -> Result<Vec<Option<ParsedSeqId>>, SeqIdError> {
    let ss = truncate_spaces(s);
    if ss.is_empty() {
        return Ok(Vec::new());
    }
    if !has_fasta_tag(ss) {
        return Ok(vec![Some(seq_id_from_text(ss, flags | NO_FASTA, warn)?)]);
    }
    let mut pieces = split_bars(ss);
    let mut ids = Vec::new();
    let mut choice = which_inverse_seq_id(pieces[0]);
    let mut tv = TypeVariant::Plain;
    if choice == SEQID_NOT_SET {
        if pieces.len() == 2 {
            choice = SEQID_GENERAL;
        }
    } else {
        tv = identify_type_variant(choice, pieces[0]);
        pieces.pop_front();
    }
    while !pieces.is_empty() {
        let step = (|| -> Result<Option<ParsedSeqId>, SeqIdError> {
            let mut id = None;
            if choice != SEQID_NOT_SET {
                let (made, next_type) = x_init(&mut pieces, choice, tv)?;
                id = made;
                choice = next_type;
            }
            if choice == SEQID_NOT_SET {
                if let Some(&front) = pieces.front() {
                    choice = which_inverse_seq_id(front);
                    if choice == SEQID_NOT_SET {
                        pieces.pop_front();
                        return Err(SeqIdError::new(
                            SEQ_ID_CPP,
                            2663,
                            "ncbi::objects::CSeq_id::ParseIDs",
                            [b"Unsupported ID type ".as_slice(), front].concat(),
                        ));
                    }
                }
            }
            if choice != SEQID_NOT_SET {
                if let Some(front) = pieces.pop_front() {
                    tv = identify_type_variant(choice, front);
                }
            }
            Ok(id)
        })();
        match step {
            Ok(id) => ids.push(id),
            Err(error) => {
                if pieces.is_empty() || flags & PARTIAL_OK == 0 {
                    return Err(error);
                }
                warn(SeqIdWarning {
                    subcode: 7,
                    message: error.what(),
                });
                while let Some(front) = pieces.front() {
                    if front.len() != 2 && front.len() != 3 {
                        pieces.pop_front();
                    } else {
                        break;
                    }
                }
            }
        }
    }
    Ok(ids)
}

#[cfg(test)]
#[path = "parse_tests.rs"]
mod tests;
