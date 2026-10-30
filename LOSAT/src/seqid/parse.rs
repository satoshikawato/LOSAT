//! The text parsers of `CSeq_id` (`c++/src/objects/seqloc/Seq_id.cpp`, `Textseq_id.cpp`,
//! `c++/src/objects/general/Object_id.cpp`) that `makeblastdb -parse_seqids` runs on a
//! definition line: `ParseIDs` (FASTA ids such as `gi|123|ref|NP_000549.1|`), the
//! constructor from a bare id (`CSeq_id(text, flags)`, with the accession guide) and the
//! typed `Set`. Each gives a `ParsedSeqId` or the `CSeqIdException` NCBI throws
//! (`SeqIdError`, whose `what` is the text NCBI prints when it posts one).

use super::accguide::{
    self, acc_type, identify_accession, parse_flags::*, SeqIdWarning, SEQID_DDBJ, SEQID_EMBL,
    SEQID_GENBANK, SEQID_GENERAL, SEQID_GI, SEQID_GIBBMT, SEQID_GIBBSQ, SEQID_GIIM, SEQID_GPIPE,
    SEQID_LOCAL, SEQID_NAMED_ANNOT_TRACK, SEQID_NOT_SET, SEQID_OTHER, SEQID_PATENT, SEQID_PDB,
    SEQID_PIR, SEQID_PRF, SEQID_SWISSPROT, SEQID_TPD, SEQID_TPE, SEQID_TPG,
};
use super::model::{
    Dbtag, GiimportId, IdPat, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId, PdbSeqId,
    TextseqId, SEQ_ID_ALTERNATIVES,
};
use crate::utils::winmask::NCBI_BUILD_ROOT;

/// A `CSeqIdException` (`eFormat`): where NCBI throws it and its message.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqIdError {
    /// The source file below `c++/src/objects/seqloc/`.
    pub file: &'static str,
    /// The line of the `NCBI_THROW`.
    pub line: u32,
    /// The function, as `NCBI_CURRENT_FUNCTION` names it.
    pub function: &'static str,
    pub message: String,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:370-384
// ```c
// string CException::ReportAll(TDiagPostFlags flags) const
// {
//     // invert the order
//     stack<const CException*> pile;
//     const CException* pex;
//     for (pex = this; pex; pex = pex->GetPredecessor()) {
//         pile.push(pex);
//     }
//     CNcbiOstrstream os;
//     os << "NCBI C++ Exception:" << '\n';
//     for (; !pile.empty(); pile.pop()) {
//         //indentation
//         os << "    ";
//         os << pile.top()->ReportThis(flags) << '\n';
//     }
// ```
// makeblastdb 2.17.0 (conda build) printed, for `>xyz|ABC|DEF`:
// `T0 "/opt/conda/conda-bld/blast_1754905306566/work/c++/src/objects/seq/../seqloc/Seq_id.cpp",
// line 2663: Error: (CSeqIdException::eFormat) ncbi::objects::CSeq_id::ParseIDs() - Unsupported
// ID type xyz`. LOSAT names the build tree of the official Linux x64 release
// (`NCBI_BUILD_ROOT`); the path is not compared (makedb-v5 DECISIONS 11).
impl SeqIdError {
    fn new(file: &'static str, line: u32, function: &'static str, message: String) -> Self {
        SeqIdError {
            file,
            line,
            function,
            message,
        }
    }

    /// `CException::what()`: the report of the exception, ending with a newline.
    pub fn what(&self) -> String {
        format!(
            "NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/objects/seq/../seqloc/{}\", line {}: Error: (CSeqIdException::eFormat) {}() - {}\n",
            self.file, self.line, self.function, self.message
        )
    }
}

const SEQ_ID_CPP: &str = "Seq_id.cpp";
const TEXTSEQ_ID_CPP: &str = "Textseq_id.cpp";

/// `isspace` in the C locale.
pub fn is_space(c: u8) -> bool {
    c == b' ' || (b'\t'..=b'\r').contains(&c)
}

/// `NStr::TruncateSpaces_Unsafe(s, NStr::eTrunc_Both)`.
pub fn truncate_spaces(s: &[u8]) -> &[u8] {
    let start = s.iter().position(|&c| !is_space(c)).unwrap_or(s.len());
    let end = s
        .iter()
        .rposition(|&c| !is_space(c))
        .map_or(start, |end| end + 1);
    &s[start..end.max(start)]
}

/// `NStr::IsBlank(s)`.
pub fn is_blank(s: &[u8]) -> bool {
    s.iter().all(|&c| is_space(c))
}

/// `NStr::StringToNonNegativeInt(str)`: the value, or -1.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:457-491
/// ```c++
/// int NStr::StringToNonNegativeInt(const CTempString str, TStringToNumFlags flags)
/// {
///     int error = 0, ret = -1;
///     size_t len = str.size();
///
///     if (!len) {
///         error = EINVAL;
///     } else {
///         size_t i = 0;
///         // skip leading '+' if any
///         if (str.data()[0] == '+'  &&  len > 1) {
///             ++i;
///         }
///         unsigned v = 0;
///         for (; i < len; ++i) {
///             unsigned d = str.data()[i] - '0';
///             if (d > 9) {
///                 error = EINVAL;
///                 break;
///             }
///             unsigned nv = v * 10 + d;
///             const unsigned kOverflowLimit = (INT_MAX - 9) / 10 + 1;
///             if (v >= kOverflowLimit) {
///                 // possible overflow
///                 if (v > kOverflowLimit || nv > INT_MAX) {
///                     error = ERANGE;
///                     break;
///                 }
///             }
///             v = nv;
///         }
///         if (!error) {
///             ret = static_cast<int>(v);
///         }
///     }
/// ```
pub fn string_to_non_negative_int(s: &[u8]) -> i32 {
    if s.is_empty() {
        return -1;
    }
    let digits = if s[0] == b'+' && s.len() > 1 {
        &s[1..]
    } else {
        s
    };
    let mut value: u32 = 0;
    for &c in digits {
        let d = u32::from(c.wrapping_sub(b'0'));
        if d > 9 {
            return -1;
        }
        let next = value.wrapping_mul(10).wrapping_add(d);
        const OVERFLOW_LIMIT: u32 = (i32::MAX as u32 - 9) / 10 + 1;
        if value >= OVERFLOW_LIMIT && (value > OVERFLOW_LIMIT || next > i32::MAX as u32) {
            return -1;
        }
        value = next;
    }
    value as i32
}

/// `NStr::StringToNumeric<TGi>(str)` (`NStr::StringToInt8` with no flags): `None` where NCBI
/// throws.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:798-871
/// ```c++
/// Int8 NStr::StringToInt8(const CTempString str, TStringToNumFlags flags, int base)
/// {
///     ...
///     // Determine sign
///     bool sign = false;
///     switch (str[pos]) {
///     case '-':
///         sign = true;
///         /*FALLTHRU*/
///     case '+':
///         pos++;
///         break;
///     ...
///     Int8 limdiv = base==10? kMax_I8 / 10: kMax_I8 / base;
///     Int8 limoff = (base==10? kMax_I8 % 10: kMax_I8 % base) + (sign ? 1 : 0);
///     ...
///     while (char ch = str[pos]) {
///         int  delta;   // corresponding numeric value of 'ch'
///         ...
///         if ( !s_IsGoodCharForRadix(ch, base, &delta) ) {
///             break;
///         }
///         // Overflow check
///         if ( n >= limdiv  &&  (n > limdiv  ||  delta > limoff) ) {
///             S2N_CONVERT_ERROR_OVERFLOW(Int8);
///         }
///         n *= base;
///         n += delta;
///         pos++;
///     }
///
///     // Last checks
///     if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
///         S2N_CONVERT_ERROR_INVAL(Int8);
///     }
///     ...
///     n = sign ? -n : n;
///     CHECK_ENDPTR(Int8);
/// ```
/// (A NUL byte ends the digits; then the end check fails.)
fn string_to_int8(s: &[u8]) -> Option<i64> {
    let (negative, digits) = match s.first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return None;
    }
    let mut n: i64 = 0;
    let limdiv = i64::MAX / 10;
    let limoff = i64::MAX % 10 + i64::from(negative);
    for &c in digits {
        let delta = i64::from(c - b'0');
        if n >= limdiv && (n > limdiv || delta > limoff) {
            return None;
        }
        n = n.wrapping_mul(10).wrapping_add(delta);
    }
    Some(if negative { n.wrapping_neg() } else { n })
}

/// `CObject_id::SetStrOrId(str)`.
///
/// NCBI reference (598d8ae6): c++/src/objects/general/Object_id.cpp:187-198
/// ```c++
/// void CObject_id::SetStrOrId(CTempString str)
/// {
///     if ( !str.empty()  &&  str[0] >= '1'  &&  str[0] <= '9' ) {
///         int id = NStr::StringToNonNegativeInt(str);
///         if ( id > 0 ) {
///             // converted to a positive integer
///             SetId(id);
///             return;
///         }
///     }
///     SetStr(str);
/// }
/// ```
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
/// ```c++
/// static const TChoiceMapEntry sc_ChoiceArray[] = {
///     { "???",          CSeq_id::e_not_set },
///     { "bbm",          CSeq_id::e_Gibbmt },
///     { "bbs",          CSeq_id::e_Gibbsq },
///     { "dbj",          CSeq_id::e_Ddbj },
///     { "emb",          CSeq_id::e_Embl },
///     { "gb",           CSeq_id::e_Genbank },
///     { "gi",           CSeq_id::e_Gi },
///     { "gibbsq",       CSeq_id::e_Gibbsq },
///     { "gim",          CSeq_id::e_Giim },
///     { "gnl",          CSeq_id::e_General },
///     { "gpp",          CSeq_id::e_Gpipe },
///     { "lcl",          CSeq_id::e_Local },
///     { "nat",          CSeq_id::e_Named_annot_track },
///     { "not_set",      CSeq_id::e_not_set },
///     { "pat",          CSeq_id::e_Patent },
///     { "pdb",          CSeq_id::e_Pdb },
///     { "pgp",          CSeq_id::e_Patent },
///     { "pir",          CSeq_id::e_Pir },
///     { "prf",          CSeq_id::e_Prf },
///     { "ref",          CSeq_id::e_Other },
///     { "sp",           CSeq_id::e_Swissprot },
///     { "tpd",          CSeq_id::e_Tpd },
///     { "tpe",          CSeq_id::e_Tpe },
///     { "tpg",          CSeq_id::e_Tpg },
///     { "tr",           CSeq_id::e_Swissprot }
/// };
/// typedef CStaticPairArrayMap<CTempString, CSeq_id::E_Choice,
///                             PNocase_Generic<CTempString> > TChoiceMap;
/// ...
/// CSeq_id::E_Choice CSeq_id::WhichInverseSeqId(const CTempString& SeqIdCode)
/// {
///     TChoiceMap::const_iterator it = sc_ChoiceMap.find(SeqIdCode);
///     if (it == sc_ChoiceMap.end()) {
///         return e_not_set;
///     } else {
///         return it->second;
///     }
/// }
/// ```
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
/// ```c++
/// static const char* const s_TextId[CSeq_id::e_MaxChoice+1] =
/// {   // FASTA_LONG formats
///     "???" , // not-set = ???
///     "lcl",  // local = lcl|integer or string
///     "bbs",  // gibbsq = bbs|integer
///     "bbm",  // gibbmt = bbm|integer
///     "gim",  // giim = gim|integer
///     "gb",   // genbank = gb|accession|locus
///     "emb",  // embl = emb|accession|locus
///     "pir",  // pir = pir|accession|name
///     "sp",   // swissprot = sp|accession|name *OR* tr|accession|name
///     "pat",  // patent = pat|country|patent number (string)|seq number (integer)
///     "ref",  // other = ref|accession|name|release - changed from oth to ref
///     "gnl",  // general = gnl|database(string)|id (string or number)
///     "gi",   // gi = gi|integer
///     "dbj",  // ddbj = dbj|accession|locus
///     "prf",  // prf = prf|accession|name
///     "pdb",  // pdb = pdb|entry name (string)|chain id (char)
///     "tpg",  // tpg = tpg|accession|name
///     "tpe",  // tpe = tpe|accession|name
///     "tpd",  // tpd = tpd|accession|name
///     "gpp",  // gpipe = gpp|accession|name
///     "nat",  // named_annot_track = nat|accession|name
///     ""  // Placeholder for end of list
/// };
/// ```
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
/// ```c++
/// inline
/// CSeq_id::ETypeVariant CSeq_id::x_IdentifyTypeVariant(CSeq_id::E_Choice type,
///                                                      const CTempString& str)
/// {
///     if (type == CSeq_id::e_Swissprot  &&  NStr::EqualNocase(str, "tr")) {
///         return eTV_tr;
///     } else if (type == CSeq_id::e_Patent  &&  NStr::EqualNocase(str, "pgp")) {
///         return eTV_pgp;
///     } else {
///         return eTV_plain;
///     }
/// }
/// ```
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
/// ```c++
/// CTextseq_id& CTextseq_id::Set
/// (const CTempString& acc_in,
///  const CTempString& name_in,
///  int                version,
///  const CTempString& release_in,
///  bool               allow_dot_version)
/// {
///     // Perform general sanity checks up front.
///     if (version < 0) {
///         NCBI_THROW(CSeqIdException, eFormat,
///                    "Unexpected negative version " + NStr::IntToString(version)
///                    + " for accession " + string(acc_in));
///     }
///
///     CTempString acc     =
///         NStr::TruncateSpaces_Unsafe(acc_in,     NStr::eTrunc_Both);
///     CTempString name    =
///         NStr::TruncateSpaces_Unsafe(name_in,    NStr::eTrunc_Both);
///     CTempString release =
///         NStr::TruncateSpaces_Unsafe(release_in, NStr::eTrunc_Both);
///     ...
///     if (acc.empty()) {
///         ResetAccession();
///     } else {
///         SIZE_TYPE idx = NPOS;
///
///         if (allow_dot_version) {
///             idx = acc.rfind('.');
///         }
///
///         if (idx == NPOS) {
///             // no version within acc
///             SetAccession(acc);
///
///             // any standalone version is ok
///             if (version > 0) {
///                 SetVersion(version);
///             } else {
///                 ResetVersion();
///             }
///         } else {
///             // accession.version
///             CTempString accession = acc.substr(0, idx);
///             CTempString acc_ver   = acc.substr(idx + 1);
///             int         ver       = NStr::StringToNonNegativeInt(acc_ver);
///
///             if (ver <= 0) {
///                 NCBI_THROW(CSeqIdException, eFormat,
///                            "Version embedded in accession " + string(acc)
///                            + " is not a positive integer");
///             } else if (version > 0  &&  ver != version) {
///                 NCBI_THROW(CSeqIdException, eFormat,
///                            "Incompatible version " + NStr::IntToString(version)
///                            + " supplied for accession " + string(acc));
///             }
///
///             SetAccession(accession);
///             SetVersion(ver);
///         }
///     }
///
///     if (name.empty()) {
///         ResetName();
///     } else {
///         SetName(name);
///     }
///
///     if (acc.empty()  &&  name.empty()) {
///         NCBI_THROW(CSeqIdException, eFormat,
///                    "Accession and name missing for Textseq-id (but got"
///                    " version " + NStr::IntToString(version) + ", release "
///                    + string(release) + ')');
///     } else if (release.empty()) {
///         ResetRelease();
///     } else {
///         SetRelease(release);
///     }
///     return *this;
/// }
/// ```
fn textseq_set(
    acc_in: &[u8],
    name_in: &[u8],
    version: i32,
    release_in: &[u8],
    allow_dot_version: bool,
) -> Result<TextseqId, SeqIdError> {
    const FUNCTION: &str = "ncbi::objects::CTextseq_id::Set";
    let lossy = |s: &[u8]| String::from_utf8_lossy(s).into_owned();
    if version < 0 {
        return Err(SeqIdError::new(
            TEXTSEQ_ID_CPP,
            67,
            FUNCTION,
            format!(
                "Unexpected negative version {version} for accession {}",
                lossy(acc_in)
            ),
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
                        format!(
                            "Version embedded in accession {} is not a positive integer",
                            lossy(acc)
                        ),
                    ));
                } else if version > 0 && ver != version {
                    return Err(SeqIdError::new(
                        TEXTSEQ_ID_CPP,
                        112,
                        FUNCTION,
                        format!(
                            "Incompatible version {version} supplied for accession {}",
                            lossy(acc)
                        ),
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
            format!(
                "Accession and name missing for Textseq-id (but got version {version}, release {})",
                lossy(release)
            ),
        ));
    } else if !release.is_empty() {
        id.release = Some(release.to_vec());
    }
    Ok(id)
}

/// `CSeq_id::Set(the_type, the_id)` of a numeric id.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2850-2884
/// ```c++
/// CSeq_id& CSeq_id::Set(E_Choice the_type, TIntId the_id)
/// {
/// // see CSeq_id::Set below, it prohibits lcl|0, but allows gi|0
///     if ((the_id < 0) || (the_type == e_Local && the_id == 0)) {
///         NCBI_THROW(CSeqIdException, eFormat,
///                    "Non-positive numeric ID " + NStr::NumericToString(the_id));
///     }
///
///     switch (the_type) {
///     case e_Local:
///         SetLocal().SetId(INT_ID_TO(CObject_id::TId, the_id));
///         break;
///     case e_Gibbsq:
///         SetGibbsq(INT_ID_TO(CSeq_id::TGibbsq, the_id));
///         break;
///     case e_Gibbmt:
///         SetGibbmt(INT_ID_TO(CSeq_id::TGibbmt, the_id));
///         break;
///     case e_Giim:
///     {
///         CGiimport_id& giim = SetGiim();
///         giim.SetId(INT_ID_TO(CGiimport_id::TId, the_id));
///         giim.ResetDb();
///         giim.ResetRelease();
///         break;
///     }
///     case e_Gi:
///         SetGi(GI_FROM(TIntId, the_id));
///         break;
///     default:
///         NCBI_THROW(CSeqIdException, eFormat,
///                    "Invalid numeric ID type" + SelectionName(the_type));
///     }
///     return *this;
/// }
/// ```
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
/// ```c++
/// CSeq_id& CSeq_id::Set(E_Choice           the_type,
///                       const CTempString& acc_in,
///                       const CTempString& name_in,
///                       int                version,
///                       const CTempString& release_in)
/// {
///     CTempString  acc       = NStr::TruncateSpaces_Unsafe(acc_in,
///                                                          NStr::eTrunc_Both);
///
///     int          the_id;
///     CTextseq_id* tsid      = 0;
///     bool         allow_dot = true;
///
///     switch (the_type) {
///     case e_not_set: // Will cause unspecified SeqId to be returned.
///         break;
///
///     case e_Local:
///         SetLocal().SetStrOrId(acc);
///         break;
///
///         // numeric IDs
///     case e_Gibbsq:
///     case e_Gibbmt:
///     case e_Giim:
///     case e_Gi:
/// #ifdef NCBI_INT8_GI
///         if ( the_type == e_Gi ) {
///             try {
///                 TGi gi = NStr::StringToNumeric<TGi>(acc);
///                 if ( gi > ZERO_GI ) {
///                     SetGi(gi);
///                     return *this;
///                 }
///             }
///             catch ( CException& /*ignored*/ ) {
///                 // will be processed by the code below
///             }
///         }
/// #endif
///         if ( (the_id = NStr::StringToNonNegativeInt (acc)) >= 0 ) {
///             return Set(the_type, the_id);
///         } else {
///             NCBI_THROW(CSeqIdException, eFormat,
///                        "Negative, excessively large, or non-numeric "
///                        + SelectionName(the_type)
///                        + " ID " + string(acc));
///         }
///         break;
///
///         // text IDs
///     case e_Genbank:    tsid = &SetGenbank();    break;
///     case e_Embl:       tsid = &SetEmbl();       break;
///     case e_Pir:        tsid = &SetPir();        allow_dot = false;  break;
///     case e_Swissprot:  tsid = &SetSwissprot();  break;
///     case e_Other:      tsid = &SetOther();      break;
///     case e_Ddbj:       tsid = &SetDdbj();       break;
///     case e_Prf:        tsid = &SetPrf();        allow_dot = false;  break;
///     case e_Tpg:        tsid = &SetTpg();        break;
///     case e_Tpe:        tsid = &SetTpe();        break;
///     case e_Tpd:        tsid = &SetTpd();        break;
///     case e_Gpipe:      tsid = &SetGpipe();      break;
///     case e_Named_annot_track:  tsid = &SetNamed_annot_track();  break;
///
///     case e_Patent:
///         {
///             CTempString      name      =
///                 NStr::TruncateSpaces_Unsafe(name_in, NStr::eTrunc_Both);
///             CTempString      release   =
///                 NStr::TruncateSpaces_Unsafe(release_in, NStr::eTrunc_Both);
///             CPatent_seq_id&  pat       = SetPatent();
///             CId_pat&         id_pat    = pat.SetCit();
///             CId_pat::C_Id&   id_pat_id = id_pat.SetId();
///             id_pat.SetCountry(acc);
///
///             if (NStr::EqualNocase(release, "pgp")) {
///                 id_pat_id.SetApp_number(name);
///             } else {
///                 id_pat_id.SetNumber(name);
///             }
///             id_pat.ResetDoc_type();
///             pat.SetSeqid(version);
///             break;
///         }
///
///     case e_General:
///         {
///             CTempString name = NStr::TruncateSpaces_Unsafe(name_in,
///                                                            NStr::eTrunc_Both);
///             CDbtag&     dbt  = SetGeneral();
///             dbt.SetDb(acc);
///             dbt.SetTag().SetStrOrId(name);
///             break;
///         }
///
///     case e_Pdb:
///         {
///             CTempString  name = NStr::TruncateSpaces_Unsafe(name_in,
///                                                             NStr::eTrunc_Both);
///             CPDB_seq_id& pdb  = SetPdb();
///             pdb.SetMol().Set(acc);
///
///             // Consult name_in in addition to name as whitespace
///             // stripping can lose relevant information here.
///             if (name_in.empty()) {
///                 pdb.ResetChain();
///             } else if (name.empty()) {
///                 pdb.SetChain(' ');
///                 name = " ";
///             } else if (name.size() == 1) {
///                 pdb.SetChain(static_cast<unsigned char>(name[0]));
///             } else {
///                 pdb.ResetChain();
///                 ERR_POST_X(16,
///                            Info << "Necessarily using backwards-incompatible"
///                            " representation for chain " << string(name)
///                            << " of PDB molecule " << acc << '.');
///             }
///             if (name.empty()) {
///                 pdb.ResetChain_id();
///             } else {
///                 pdb.SetChain_id(name);
///             }
///             pdb.ResetRel();
///             break;
///         }
///
///     default:
///         NCBI_THROW(CSeqIdException, eFormat,
///                    "Unsupported Seq-id type " + SelectionName(the_type));
///     }
///
///     if (tsid) {
///         // CTextseq_id::Set will take care of truncating any spaces.
///         tsid->Set(acc, name_in, version, release_in, allow_dot);
///     }
///
///     return *this;
/// }
/// ```
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
            if choice == SEQID_GI {
                if let Some(gi) = string_to_int8(acc) {
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
                    format!(
                        "Negative, excessively large, or non-numeric {} ID {}",
                        selection_name(choice),
                        String::from_utf8_lossy(acc)
                    ),
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
/// ```c++
/// static void s_SplitVersion(const CTempString& acc_in, CTempString& acc,
///                            int& ver)
/// {
///     CTempString verstr;
///     NStr::SplitInTwo(acc_in, ".", acc, verstr);
///     if (verstr.empty()) {
///         ver = 0;
///     } else {
///         ver = NStr::StringToNonNegativeInt(verstr);
///         if (ver <= 0) {
///             NCBI_THROW(CSeqIdException, eFormat,
///                        "Version embedded in accession " + string(acc_in)
///                        + " is not a positive integer");
///         }
///     }
/// }
/// ```
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
            format!(
                "Version embedded in accession {} is not a positive integer",
                String::from_utf8_lossy(acc_in)
            ),
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
/// ```c++
/// bool CSeq_id::IsValidLocalID(const CTempString& s)
/// {
///     return (fNoError == CheckLocalID(s));
/// }
///
///
/// CSeq_id::TErrorFlags
/// CSeq_id::CheckLocalID(const CTempString& s)
/// {
///     if (NStr::IsBlank(s)) {
///         return fEmptyId;
///     }
///
///     TErrorFlags error_flags = fNoError;
///     if (s.length() > kMaxLocalIDLength) {
///         error_flags |= fExceedsMaxLength;
///     }
///
///     static const char* kIllegal = " >[]|\"";
///     CSeq_id_find_pred pred; pred.kSymbols = kIllegal;
///     if (find_if(s.begin(), s.end(), pred) != s.end()) {
///         error_flags |= fInvalidChar;
///     }
///     return error_flags;
/// }
/// ```
pub fn is_valid_local_id(s: &[u8]) -> bool {
    !is_blank(s) && s.len() <= MAX_LOCAL_ID_LENGTH && !has_invalid_local_char(s)
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:77-86
// ```c++
// struct CSeq_id_find_pred
// {
//     const char* kSymbols;
//     bool operator()(const char ch)
//     {
//         if (ch<32 || ch>127) // non-ASCII is not allowed
//             return true;
//         return strchr(kSymbols, ch) != 0;
//     }
// };
// ```
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
/// ```c++
///     const static size_t kMaxLocalIDLength    = 50;
///     const static size_t kMaxGeneralDBLength  = 20;
///     const static size_t kMaxGeneralTagLength = 50;
///     const static size_t kMaxAccessionLength  = 30;
/// ```
pub const MAX_LOCAL_ID_LENGTH: usize = 50;
pub const MAX_GENERAL_TAG_LENGTH: usize = 50;
pub const MAX_ACCESSION_LENGTH: usize = 30;

/// `kSupportedRawDbtags`: the databases of a bare `db:tag` that is a general id.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:96-124
/// ```c++
/// static const char* sc_SupportedRawDbtags[] = {
///     "ATGC",
///     "BCMHGSC",
///     "BERKELEY",
///     "CELERA",
///     "GSDB",
///     "HOOD",
///     "LANLCHGS",
///     "LRG",
///     "MIPS",
///     "NCBI_EXT_ACC",
///     "NCBI_GENOMES",
///     "NCBI_MITO",
///     "PGEC",
///     "PID",
///     "SGD",
///     "SHGC",
///     "SRA",
///     "TIGR",
///     "UOKNOR",
///     "UWGC",
///     "WASHU",
///     "WIBR",
///     "WUGSC",
///     "dbGSS",
///     "dbSTS"
/// };
/// DEFINE_STATIC_ARRAY_MAP_WITH_COPY(CStaticArraySet<string>, kSupportedRawDbtags,
///                                   sc_SupportedRawDbtags);
/// ```
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
/// ```c++
/// static CSeq_id::E_Choice s_CheckForFastaTag(const CTempString& s)
/// {
///     // > rather than >= because there should be content after the bar.
///     if (s.size() > 3  &&  s[2] == '|') {
///         return CSeq_id::WhichInverseSeqId(s.substr(0, 2));
///     } else if (s.size() > 4  &&  s[3] == '|') {
///         return CSeq_id::WhichInverseSeqId(s.substr(0, 3));
///     } else {
///         return CSeq_id::e_not_set;
///     }
/// }
/// ```
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
/// ```c++
/// CSeq_id& CSeq_id::Set(const CTempString& the_id_in, TParseFlags flags)
/// {
///     CTempString the_id = NStr::TruncateSpaces_Unsafe(the_id_in,
///                                                      NStr::eTrunc_Both);
///     E_Choice    type   = e_not_set;
///
///     if ((flags & fParse_NoFASTA) == 0) {
///         type = s_CheckForFastaTag(the_id);
///     }
///     if (type == e_not_set) {
///         if (the_id.empty()) {
///             NCBI_THROW(CSeqIdException, eFormat,
///                        "Empty bare accession supplied");
///         }
///         // If no (attempt at a) valid tag, tries to interpret the string
///         // as a pure accession.
///         if ((flags & fParse_AnyRaw) != 0) {
///             type = GetAccType(IdentifyAccession(the_id,
///                                                 flags | fParse_FallbackOK));
///         }
///         switch (type) {
///         case e_Gi:
///             return Set(type, the_id);
///         case e_not_set:
///         {
///             // Check for general IDs, albeit only with well-known
///             // database names like SRA.
///             SIZE_TYPE colon_pos = the_id.find(':');
///             if (colon_pos != NPOS) {
///                 string db = the_id.substr(0, colon_pos);
///                 NStr::ToUpper(db);
///                 // const auto& whitelist = (*s_Guide)->general;
///                 const auto& whitelist = kSupportedRawDbtags;
///                 if (whitelist.find(db) != whitelist.end()) {
///                     // Reextract prefix to preserve case.
///                     return Set(e_General, the_id.substr(0, colon_pos),
///                                the_id.substr(colon_pos + 1));
///                 }
///             }
///             if ((flags & fParse_ValidLocal) != 0
///                 &&  ((flags & fParse_AnyLocal) == fParse_AnyLocal
///                      ||  IsValidLocalID(the_id))) {
///                 return Set(e_Local, the_id);
///             } else {
///                 NCBI_THROW(CSeqIdException, eFormat,
///                            "Malformatted ID " + string(the_id));
///             }
///         }
///         case e_Prf:
///             // technically a name/locus, not an accession!
///             return Set(type, kEmptyStr, the_id);
///         case e_Pdb:
///         {
///             string mol(the_id, 0, 4), chain;
///             // NStr::SplitInTwo(the_id, "|", mol, chain);
///             if (the_id.size() > 5) {
///                 chain = the_id.substr(5);
///             } else if (the_id.size() == 5  &&  the_id[4] != '|') {
///                 chain = the_id[4];
///             }
///             return Set(type, mol, chain);
///         }
///         default:
///         {
///             CTempString acc;
///             int         ver;
///             s_SplitVersion(the_id, acc, ver);
///             return Set(type, acc, kEmptyStr, ver);
///         }
///         }
///     } else {
///         list<CTempString> fasta_pieces;
///         NStr::Split(the_id, "|", fasta_pieces);
///         ETypeVariant tv = x_IdentifyTypeVariant(type, fasta_pieces.front());
///         fasta_pieces.pop_front();
///         x_Init(fasta_pieces, type, tv);
///         if ( !fasta_pieces.empty() ) {
///             // tolerate trailing parts if they're all empty.
///             ITERATE(list<CTempString>, it, fasta_pieces) {
///                 if ( !it->empty() ) {
///                     if ((flags & fParse_PartialOK) != 0) {
///                         ERR_POST_X(10, Warning << "Ignoring extra parts"
///                                    " (synonyms?) in FASTA-style ID "
///                                    << the_id);
///                     } else {
///                         NCBI_THROW(CSeqIdException, eFormat,
///                                    "FASTA-style ID " + string(the_id)
///                                    + " has too many parts.");
///                     }
///                 }
///             }
///         }
///         return *this;
///     }
/// }
/// ```
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
    let lossy = |s: &[u8]| String::from_utf8_lossy(s).into_owned();
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
                        format!("Malformatted ID {}", lossy(the_id)),
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
                    message: format!(
                        "Ignoring extra parts (synonyms?) in FASTA-style ID {}",
                        lossy(the_id)
                    ),
                });
            } else {
                return Err(SeqIdError::new(
                    SEQ_ID_CPP,
                    2542,
                    FUNCTION,
                    format!("FASTA-style ID {} has too many parts.", lossy(the_id)),
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
/// ```c++
/// CSeq_id::E_Choice CSeq_id::x_Init(list<CTempString>& fasta_pieces,
///                                   E_Choice type, ETypeVariant tv)
/// {
///     _ASSERT(!fasta_pieces.empty());
///     _ASSERT(type != e_not_set);
///
///     vector<CTempString> fields(3);
///     SIZE_TYPE   min_fields, max_fields;
///     E_Choice    next_type = e_not_set;
///     switch (type) {
///     case e_Local:
///     case e_Gibbsq:
///     case e_Gibbmt:
///     case e_Giim:
///     case e_Gi:
///         min_fields = max_fields = 1;
///         break;
///     case e_Patent:
///         min_fields = max_fields = 3;
///         break;
///     case e_General:
///         min_fields = max_fields = 2;
///         break;
///     ...
///     default: // text seqid: accession and optional name
///         min_fields = 1;
///         max_fields = 2;
///         break;
///     }
///
///     for (SIZE_TYPE i = 0;  i < max_fields;  ++i) {
///         if (fasta_pieces.empty()) {
///             if (i >= min_fields) {
///                 break;
///             } else {
///                 NCBI_THROW(CSeqIdException, eFormat,
///                            "Not enough fields for ID of type "
///                            + string(s_TextId[type]));
///             }
///         } else {
///             if (i >= min_fields  &&  fasta_pieces.size() > 1
///                 &&  (fasta_pieces.front().size() == 2
///                      ||  fasta_pieces.front().size() == 3)
///                 &&  ((next_type = WhichInverseSeqId(fasta_pieces.front()))
///                      != e_not_set)) {
///                 // Likely mid-string optional-field omission; look ahead
///                 // more to see whether the piece works better as a field
///                 // or a tag.
///                 list<CTempString>::iterator it = fasta_pieces.begin();
///                 ++it;
///                 _ASSERT(it != fasta_pieces.end());
///                 E_Choice next_type_2;
///                 if ((it->size() == 2  ||  it->size() == 3)
///                     &&  (next_type_2 = WhichInverseSeqId(*it)) != e_not_set) {
///                     next_type = next_type_2;
///                 } else if (it->find_first_not_of(" \t\n") == NPOS
///                            &&  ++it == fasta_pieces.end()) {
///                     next_type = e_not_set;
///                 } else {
///                     break;
///                 }
///             }
///             fields[i] = fasta_pieces.front();
///             fasta_pieces.pop_front();
///         }
///     }
///
///     // Special case -- dbSNP IDs have historically contained internal
///     // vertical bars, so we have to parse them greedily.
///     string snp_name; // must survive until the end of the function
///     if (type == e_General  &&  NStr::EqualNocase(fields[0], "dbSNP")
///         &&  !fasta_pieces.empty() ) {
///         snp_name = string(fields[1]) + '|' + NStr::Join(fasta_pieces, "|");
///         fields[1] = snp_name;
///         fasta_pieces.clear();
///     }
///
///     // Clear out extra empty pieces
///     while ( !fasta_pieces.empty()  &&  fasta_pieces.front().empty() ) {
///         fasta_pieces.pop_front();
///     }
///
///     int ver = 0;
///     switch (type) {
///     case e_Swissprot:
///         if (tv == eTV_tr) {
///             fields[2] = "unreviewed";
///         } else {
///             fields[2] = "reviewed";
///         }
///         break;
///
///     case e_Patent:
///         // "version" actually sequence number within patent, but whatever...
///         ver = NStr::StringToNonNegativeInt(fields[2]);
///         if (ver < 0) {
///             NCBI_THROW(CSeqIdException, eFormat,
///                        "Bad sequence number " + string(fields[2]) + " for "
///                        + string(fields[0]) + " patent " + string(fields[1]));
///         }
///         // to distinguish applications from granted patents; the numeric
///         // content has already made its way into ver.
///         fields[2] = "pat";
///         break;
///
///     case e_Pdb:
///         if (fields[0].size() < 4
///             ||  (fields[0].size() > 5
///                  &&  ( !fields[1].empty()
///                        ||  strchr("|-_", fields[0][4]) == NULL))) {
///             NCBI_THROW(CSeqIdException, eFormat,
///                        "Malformatted PDB ID " + string(fields[0]));
///         }
///         if (fields[0].size() > 4  &&  fields[1].empty()) { // misdelimited
///             if (fields[0].size() > 5) {
///                 fields[1] = fields[0].substr(5);
///             } else {
///                 _ASSERT(fields[0][4] != '|');
///                 fields[1] = fields[0].substr(4);
///             }
///             fields[0] = fields[0].substr(0, 4);
///         }
///         break;
///
///     default:
///         break; // avoid compiler warnings
///     }
///
///     Set(type, fields[0] /* acc */, fields[1] /* name */, ver,
///         fields[2] /* rel */);
///
///     return next_type;
/// }
/// ```
/// (`strchr` also finds the NUL of a `fields[0]` whose fifth byte is NUL.)
fn x_init(
    pieces: &mut std::collections::VecDeque<&[u8]>,
    choice: u32,
    tv: TypeVariant,
) -> Result<(Option<ParsedSeqId>, u32), SeqIdError> {
    const FUNCTION: &str = "ncbi::objects::CSeq_id::x_Init";
    let lossy = |s: &[u8]| String::from_utf8_lossy(s).into_owned();
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
                    format!(
                        "Bad sequence number {} for {} patent {}",
                        lossy(&fields[2]),
                        lossy(&fields[0]),
                        lossy(&fields[1])
                    ),
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
                    format!("Malformatted PDB ID {}", lossy(mol)),
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
/// ```c++
/// SIZE_TYPE CSeq_id::ParseIDs(CBioseq::TId& ids, const CTempString& s,
///                             TParseFlags flags)
/// {
///     CTempString ss = NStr::TruncateSpaces_Unsafe(s, NStr::eTrunc_Both);
///     if (ss.empty()) {
///         return 0;
///     }
///
///     // first simple check to make it faster
///     if (!s_HasFastaTag(ss)) {
///         CRef<CSeq_id> id(new CSeq_id(ss, flags | fParse_NoFASTA));
///         ids.push_back(id);
///         return 1;
///     }
///
///     SIZE_TYPE count = 0;
///     list<CTempString> fasta_pieces;
///     NStr::Split(ss, "|", fasta_pieces);
///     _ASSERT(fasta_pieces.size() > 0);
///     if (fasta_pieces.size() == 1)
///     {
///         CRef<CSeq_id> id(new CSeq_id(ss, flags | fParse_NoFASTA));
///         ids.push_back(id);
///         count = 1;
///     }
///     else
///     {
///         E_Choice     type = WhichInverseSeqId(fasta_pieces.front());
///         ETypeVariant tv;
///         if (type == e_not_set) {
///             if (fasta_pieces.size() == 2) {
///                 // unknown database are reported as 'general'
///                 type = e_General;
///             }
///             tv   = eTV_plain;
///         } else {
///             tv = x_IdentifyTypeVariant(type, fasta_pieces.front());
///             fasta_pieces.pop_front();
///         }
///         while ( !fasta_pieces.empty() ) {
///             try {
///                 CRef<CSeq_id> id(new CSeq_id);
///                 if (type != e_not_set) {
///                     type = id->x_Init(fasta_pieces, type, tv);
///                 }
///                 if (type == e_not_set  &&  !fasta_pieces.empty() ) {
///                     type = WhichInverseSeqId(fasta_pieces.front());
///                     if (type == e_not_set) {
///                         CTempString typestr = fasta_pieces.front();
///                         fasta_pieces.pop_front();
///                         NCBI_THROW(CSeqIdException, eFormat,
///                                    "Unsupported ID type " + typestr);
///                     }
///                 }
///                 if (type != e_not_set) {
///                     _ASSERT( !fasta_pieces.empty() );
///                     tv = x_IdentifyTypeVariant(type, fasta_pieces.front());
///                     fasta_pieces.pop_front();
///                 }
///                 ids.push_back(id);
///                 ++count;
///             } catch (std::exception& e) {
///                 if (fasta_pieces.empty()) {
///                     throw;
///                 }
///                 if ((flags & fParse_PartialOK) != 0) {
///                     ERR_POST_X(7, Warning << e.what());
///                     do {
///                         auto l = fasta_pieces.front().size();
///                         if (l != 2  &&  l != 3) {
///                             fasta_pieces.pop_front();
///                         } else {
///                             break;
///                         }
///                     } while ( !fasta_pieces.empty() );
///                 } else {
///                     throw;
///                 }
///             }
///         }
///     }
///     return count;
/// }
/// ```
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
                            format!("Unsupported ID type {}", String::from_utf8_lossy(front)),
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

/// The accession guide, for the callers that classify ids (`IdentifyAccession`).
pub use accguide::identify_seq_id;

#[cfg(test)]
#[path = "parse_tests.rs"]
mod tests;
