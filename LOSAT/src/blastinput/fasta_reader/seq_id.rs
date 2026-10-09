//! The first line of a record that `CBlastInputReader` tries as a Seq-id.
//!
//! NCBI reads such a line as the ID of a sequence to fetch with its data loaders (GenBank
//! over the network, or a BLAST database), which LOSAT does not do; it is an explicit
//! rejection. A line that NCBI's `CSeq_id` parser rejects as a malformatted ID is read as
//! FASTA data, as NCBI does.
//!
//! The parse is `CSeq_id::Set` (Seq_id.cpp:2457-2551) with `fParse_AnyRaw |
//! fParse_ValidLocal`, then the retry with `fParse_AnyRaw` when the result is a local ID.
//! Every branch is ported except the accession guide (`SAccGuide::Find` over the 1458
//! prefix rules and the special ranges of `accguide2.inc`): a line in one of the guide's
//! formats is rejected without the lookup (`SeqIdLine::GuideFormat`), which also rejects
//! some lines that the guide does not know and NCBI reads as FASTA (`ZZ123456`).

use super::ReaderConfig;

/// What `CSeq_id(line, fParse_AnyRaw | fParse_ValidLocal)` and its retry make of a line.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum SeqIdLine {
    /// `Malformatted ID`: `CFastaReader` reads the line as data.
    Fasta,
    /// A FASTA tag (`gb|`, `lcl|`, ...): a Seq-id, or a `CSeqIdException` that is not
    /// `Malformatted ID` (rethrown).
    FastaTag,
    /// A GI (`eAcc_gi`).
    Gi,
    /// A PDB ID (`eAcc_pdb`).
    Pdb,
    /// A Swiss-Prot accession (`eAcc_swissprot`).
    Swissprot,
    /// A PRF name (`eAcc_prf`).
    Prf,
    /// Letters (or `_`) and digits in a format that the accession guide lists: an
    /// accession unless the guide's rules for that prefix say otherwise.
    GuideFormat,
    /// `<db>:<tag>` with a database of the raw general-ID whitelist.
    General,
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:526-568
// ```c++
// static const TChoiceMapEntry sc_ChoiceArray[] = {
//     { "???",          CSeq_id::e_not_set },
//     { "bbm",          CSeq_id::e_Gibbmt },
//     { "bbs",          CSeq_id::e_Gibbsq },
//     { "dbj",          CSeq_id::e_Ddbj },
//     { "emb",          CSeq_id::e_Embl },
//     { "gb",           CSeq_id::e_Genbank },
//     { "gi",           CSeq_id::e_Gi },
//     { "gibbsq",       CSeq_id::e_Gibbsq },
//     { "gim",          CSeq_id::e_Giim },
//     { "gnl",          CSeq_id::e_General },
//     { "gpp",          CSeq_id::e_Gpipe },
//     { "lcl",          CSeq_id::e_Local },
//     { "nat",          CSeq_id::e_Named_annot_track },
//     { "not_set",      CSeq_id::e_not_set },
//     { "pat",          CSeq_id::e_Patent },
//     { "pdb",          CSeq_id::e_Pdb },
//     { "pgp",          CSeq_id::e_Patent },
//     { "pir",          CSeq_id::e_Pir },
//     { "prf",          CSeq_id::e_Prf },
//     { "ref",          CSeq_id::e_Other },
//     { "sp",           CSeq_id::e_Swissprot },
//     { "tpd",          CSeq_id::e_Tpd },
//     { "tpe",          CSeq_id::e_Tpe },
//     { "tpg",          CSeq_id::e_Tpg },
//     { "tr",           CSeq_id::e_Swissprot }
// };
// typedef CStaticPairArrayMap<CTempString, CSeq_id::E_Choice,
//                             PNocase_Generic<CTempString> > TChoiceMap;
// ```
// (commented-out aliases omitted). The tags that give a type other than `e_not_set` and
// that a 2- or 3-byte prefix can match; `gibbsq` is six bytes.
const FASTA_TAGS: [&[u8]; 22] = [
    b"bbm", b"bbs", b"dbj", b"emb", b"gb", b"gi", b"gim", b"gnl", b"gpp", b"lcl", b"nat", b"pat",
    b"pdb", b"pgp", b"pir", b"prf", b"ref", b"sp", b"tpd", b"tpe", b"tpg", b"tr",
];

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:96-124
// ```c++
// static const char* sc_SupportedRawDbtags[] = {
//     "ATGC",
//     "BCMHGSC",
//     "BERKELEY",
//     "CELERA",
//     "GSDB",
//     "HOOD",
//     "LANLCHGS",
//     "LRG",
//     "MIPS",
//     "NCBI_EXT_ACC",
//     "NCBI_GENOMES",
//     "NCBI_MITO",
//     "PGEC",
//     "PID",
//     "SGD",
//     "SHGC",
//     "SRA",
//     "TIGR",
//     "UOKNOR",
//     "UWGC",
//     "WASHU",
//     "WIBR",
//     "WUGSC",
//     "dbGSS",
//     "dbSTS"
// };
// DEFINE_STATIC_ARRAY_MAP_WITH_COPY(CStaticArraySet<string>, kSupportedRawDbtags,
//                                   sc_SupportedRawDbtags);
// ```
// A `CStaticArraySet<string>` compares case-sensitively, and the key is upper-cased
// before the lookup, so `dbGSS` and `dbSTS` never match.
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

// NCBI reference (598d8ae6): c++/src/objects/seqloc/accguide2.inc:36-49
// ```c++
// static const char* const kBuiltInGuide[] = {
//     "# $Id: accguide2.inc 696906 2025-04-28 18:48:20Z ivanov $",
//     "version  2 # of file format",
//     "",
//     "# three-letter-prefix protein accessions (traditionally with five digits)",
//     "3+5  AAE  gb_patent_prot",
//     "3+5  ??_  unknown # Longer variants (6-9 digits as of June 2018) are RefSeq",
//     "3+5  A??  gb_prot",
//     "3+7  A??  gb_prot",
//     "3+9  A??  gb_prot",
//     "3+11 A??  gb_prot",
//     "3+5  B??  ddbj_prot",
//     "3+7  B??  ddbj_prot",
//     "3+9  B??  ddbj_prot",
// ```
// The formats (`<letters>+<digits>`) of the 1458 prefix rules of `kBuiltInGuide`
// (lines 36-25779; for example 1+5 from line 319, 2+6 from 348, 3+6 from 277, 4+5 from
// 290); `SAccGuide::Find` returns `eAcc_unknown` for any other format
// (Seq_id.cpp:1374-1378). Every rule's prefix is made of `A`-`Z`, `_` and the wildcard
// `?`.
const GUIDE_FORMATS: [(usize, usize); 26] = [
    (1, 5),
    (2, 6),
    (2, 8),
    (2, 10),
    (3, 5),
    (3, 6),
    (3, 7),
    (3, 8),
    (3, 9),
    (3, 11),
    (4, 5),
    (4, 6),
    (4, 8),
    (4, 9),
    (4, 10),
    (5, 6),
    (5, 7),
    (6, 9),
    (6, 10),
    (6, 11),
    (7, 8),
    (7, 9),
    (7, 10),
    (9, 9),
    (9, 10),
    (9, 11),
];

/// `isalnum`, `isalpha` and `isdigit` in the C locale.
fn alnum(byte: u8) -> bool {
    byte.is_ascii_alphanumeric()
}

fn alpha(byte: u8) -> bool {
    byte.is_ascii_alphabetic()
}

fn digit(byte: u8) -> bool {
    byte.is_ascii_digit()
}

/// The type of a FASTA tag at the start of the line, if it is one of `FASTA_TAGS`.
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
fn has_fasta_tag(s: &[u8]) -> bool {
    let tag = if s.len() > 3 && s[2] == b'|' {
        &s[..2]
    } else if s.len() > 4 && s[3] == b'|' {
        &s[..3]
    } else {
        return false;
    };
    FASTA_TAGS
        .iter()
        .any(|known| known.eq_ignore_ascii_case(tag))
}

/// The accession type of a line without a FASTA tag (`None` for `eAcc_unknown` and for
/// the types that `GetAccType` maps to `e_not_set`; `GuideFormat` for every guide format).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1634-1663
/// ```c++
/// CSeq_id::EAccessionInfo CSeq_id::IdentifyAccession(const CTempString& acc,
///                                                    TParseFlags flags)
/// {
///     SIZE_TYPE main_size = acc.find('.');
///     bool has_version = true;
///     if (main_size == NPOS) {
///         has_version = false;
///         main_size = acc.size();
///     } else if (main_size >= acc.size() - 1
///                ||  acc.find_first_not_of(kDigits, main_size + 1) != NPOS) {
///         return eAcc_unknown; // non-numeric "version"
///     }
///
///     static const SIZE_TYPE kMainAccBufSize = 32;
///     if (main_size <= kMainAccBufSize) {
///         const unsigned char* ucdata = (const unsigned char*)acc.data();
///         char main_acc_buf[kMainAccBufSize];
///         for (SIZE_TYPE i = 0;  i < main_size;  ++i) {
///             main_acc_buf[i] = toupper(ucdata[i]);
///         }
///         CTempString main_acc(main_acc_buf, main_size);
///         return x_IdentifyAccession(main_acc, flags, has_version);
///     } else {
///         // Unlikely to prove recognizable (far too long for any standard
///         // format as of January 2016), but try anyway.
///         string main_acc(acc, 0, main_size);
///         NStr::ToUpper(main_acc);
///         return x_IdentifyAccession(main_acc, flags, has_version);
///     }
/// }
/// ```
fn identify_accession(acc: &[u8]) -> Option<SeqIdLine> {
    let (main, has_version) = match acc.iter().position(|&byte| byte == b'.') {
        None => (acc, false),
        Some(dot) => {
            if dot + 1 >= acc.len() || !acc[dot + 1..].iter().all(|&byte| digit(byte)) {
                return None;
            }
            (&acc[..dot], true)
        }
    };
    x_identify_accession(&main.to_ascii_uppercase(), has_version)
}

/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1665-1677
/// ```c++
/// CSeq_id::EAccessionInfo
/// CSeq_id::x_IdentifyAccession(const CTempString& main_acc, TParseFlags flags,
///                              bool has_version)
/// {
///     SIZE_TYPE digit_pos = main_acc.find_first_of(kDigits),
///         main_size = main_acc.size();
///     char flag_char = '\0';
///     if (digit_pos == NPOS) {
///         return eAcc_unknown;
///     } else {
///         SIZE_TYPE non_dig_pos = main_acc.find_first_not_of(kDigits, digit_pos);
///         const unsigned char* ucdata = (const unsigned char*)main_acc.data();
///         if (non_dig_pos != NPOS  &&  (flags & fParse_RawText) != 0) {
/// ```
/// The flags are `fParse_AnyRaw | fParse_ValidLocal | fParse_FallbackOK`
/// (`fParse_RawText` and `fParse_RawGI` set).
fn x_identify_accession(main: &[u8], has_version: bool) -> Option<SeqIdLine> {
    let digit_pos = main.iter().position(|&byte| digit(byte))?;
    let main_size = main.len();
    let mut flag_char = false;
    let non_dig_pos = main[digit_pos..]
        .iter()
        .position(|&byte| !digit(byte))
        .map(|offset| digit_pos + offset);
    if let Some(non_dig_pos) = non_dig_pos {
        // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1678-1711
        // ```c++
        //             if ( !has_version  &&  digit_pos == 0  &&  main_size >= 4
        //                 &&  non_dig_pos < 5  &&  isalnum(ucdata[1])
        //                 &&  isalnum(ucdata[2])  &&  isalnum(ucdata[3])) {
        //                 // Possible PDB (always unversioned); examine further
        //                 // to avoid false positives.
        //                 if (main_size > 4  &&  main_size <= 17
        //                     &&  strchr("|-_", main_acc[4])
        //                     &&  (main_size <= 6  ||  isalnum(ucdata[5]))) {
        //                     // Conventionally delimited
        //                     return eAcc_pdb;
        //                 } else switch (main_size) {
        //                 ...
        //                 case 4:
        //                     return eAcc_pdb;
        //                 }
        //             }
        // ```
        // (NCBI has the cases 7, 6 and 5 inside a `/* */` comment; the elided lines are that comment.)
        // `strchr` also finds the terminating NUL, so a NUL byte at index 4 counts as a
        // delimiter.
        if !has_version
            && digit_pos == 0
            && main_size >= 4
            && non_dig_pos < 5
            && alnum(main[1])
            && alnum(main[2])
            && alnum(main[3])
        {
            if main_size > 4
                && main_size <= 17
                && matches!(main[4], b'|' | b'-' | b'_' | 0)
                && (main_size <= 6 || alnum(main[5]))
            {
                return Some(SeqIdLine::Pdb);
            } else if main_size == 4 {
                return Some(SeqIdLine::Pdb);
            }
        }
        // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1713-1750
        // ```c++
        //             if (digit_pos == 1  &&  main_size == 6
        //                 &&  (main_acc[0] == 'O'  ||  main_acc[0] == 'P'
        //                      ||  main_acc[0] == 'Q' ||  isalpha(ucdata[2]))
        //                 &&  isdigit(ucdata[1])  &&  isalnum(ucdata[2])
        //                 &&  isalnum(ucdata[3])  &&  isalnum(ucdata[4])
        //                 &&  isdigit(ucdata[5])) {
        //                 return eAcc_swissprot;
        //             } else if (digit_pos == 1  &&  main_size == 10
        //                        &&  main_acc[0] != 'O'  &&  main_acc[0] != 'P'
        //                        &&  main_acc[0] != 'Q'
        //                        &&  isalpha(ucdata[2])  &&  isalnum(ucdata[3])
        //                        &&  isalnum(ucdata[4])  &&  isdigit(ucdata[5])
        //                        &&  isalpha(ucdata[6])  &&  isalnum(ucdata[7])
        //                        &&  isalnum(ucdata[8])  &&  isdigit(ucdata[9])) {
        //                 return eAcc_swissprot;
        //             } else if ( !has_version  &&  digit_pos == 0
        //                        &&  (non_dig_pos == 6  ||  non_dig_pos == 7)
        //                        &&  (main_size == non_dig_pos + 1
        //                             ||  main_acc[non_dig_pos + 1] == ':'
        //                             ||  (isalpha(ucdata[non_dig_pos + 1])
        //                                  &&  (main_size == non_dig_pos + 2
        //                                       ||  main_acc[non_dig_pos + 2] == ':')))) {
        //                 // A formal spec appears to be elusive, but all examples in ID
        //                 // contain six or seven digits followed by one or two letters,
        //                 // followed in some rare cases by a tag such as :PDB=...
        //                 return eAcc_prf;
        //             } else if (digit_pos >= 4  &&  non_dig_pos == digit_pos + 2
        //                        &&  main_size - non_dig_pos >= 6  &&  main_acc[3] != '_'
        //                        &&  (main_acc[non_dig_pos] == 'S'
        //                             ||  main_acc[non_dig_pos] == 'P')
        //                        &&  (main_acc.find_first_not_of
        //                             (kDigits, non_dig_pos + 1) == NPOS)) {
        //                 flag_char = main_acc[non_dig_pos];
        //             } else {
        //                 return eAcc_unknown;
        //             }
        //         }
        //     }
        // ```
        let oqp = matches!(main[0], b'O' | b'P' | b'Q');
        if digit_pos == 1
            && main_size == 6
            && (oqp || alpha(main[2]))
            && digit(main[1])
            && alnum(main[2])
            && alnum(main[3])
            && alnum(main[4])
            && digit(main[5])
        {
            return Some(SeqIdLine::Swissprot);
        } else if digit_pos == 1
            && main_size == 10
            && !oqp
            && alpha(main[2])
            && alnum(main[3])
            && alnum(main[4])
            && digit(main[5])
            && alpha(main[6])
            && alnum(main[7])
            && alnum(main[8])
            && digit(main[9])
        {
            return Some(SeqIdLine::Swissprot);
        } else if !has_version
            && digit_pos == 0
            && (non_dig_pos == 6 || non_dig_pos == 7)
            && (main_size == non_dig_pos + 1
                || main[non_dig_pos + 1] == b':'
                || (alpha(main[non_dig_pos + 1])
                    && (main_size == non_dig_pos + 2 || main[non_dig_pos + 2] == b':')))
        {
            return Some(SeqIdLine::Prf);
        } else if digit_pos >= 4
            && non_dig_pos == digit_pos + 2
            && main_size - non_dig_pos >= 6
            && main[3] != b'_'
            && matches!(main[non_dig_pos], b'S' | b'P')
            && main[non_dig_pos + 1..].iter().all(|&byte| digit(byte))
        {
            flag_char = true;
        } else {
            return None;
        }
    }
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1752-1769
    // ```c++
    //     if (digit_pos == 0) {
    //         if ((flags & fParse_RawGI) != 0  &&  !has_version
    //             &&  main_acc[0] != '0'
    //             &&  main_acc.find_first_not_of(kDigits) == NPOS) {
    //             return eAcc_gi; // just digits
    //         } else {
    //             return eAcc_unknown; // PDB already handled
    //         }
    //     } else if ((flags & fParse_RawText) == 0) {
    //         return eAcc_unknown;
    //     }
    //
    //     SIZE_TYPE flag_len = (flag_char == '\0') ? 0 : 1;
    //     SIZE_TYPE digit_count = main_size - digit_pos - flag_len;
    //     auto& guide = *s_Guide;
    //     const EAccessionInfo& found_ai
    //         = guide->Find(SAccGuide::s_Key(digit_pos, digit_count), main_acc);
    //     EAccessionInfo ai = found_ai;
    // ```
    if digit_pos == 0 {
        return (!has_version && main[0] != b'0' && main.iter().all(|&byte| digit(byte)))
            .then_some(SeqIdLine::Gi);
    }
    let digit_count = main_size - digit_pos - usize::from(flag_char);
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1370-1398
    // ```c++
    // const SAccGuide::TAccInfo& SAccGuide::Find(TFormatCode fmt,
    //                                            const CTempString& acc_or_pfx,
    //                                            string* key_used) const
    // {
    //     static const TAccInfo kUnknown = CSeq_id::eAcc_unknown;
    //     TMainMap::const_iterator it = rules.find(fmt);
    //     if (it == rules.end()) {
    //         return kUnknown;
    //     }
    //
    //     const SSubMap&            submap = it->second;
    //     const TAccInfo*           result = &kUnknown;
    //     CTempString               pfx     (acc_or_pfx, 0, fmt >> 16);
    //     TPrefixes::const_iterator pit    = submap.prefixes.find(pfx);
    //     if (pit != submap.prefixes.end()) {
    //         result = &pit->second;
    //     } else {
    //         ITERATE (TPairs, wit, submap.wildcards) {
    //             if (NStr::MatchesMask(pfx, wit->first)) {
    //                 bool bad_match = false; // Limit ? to matching letters
    //                 SIZE_TYPE pos = wit->first.find('?');
    //                 while (pos != NPOS) {
    //                     if ( !isalnum(pfx[pos])  &&  pfx[pos] != '?' ) {
    //                         bad_match = true;
    //                         break;
    //                     } else {
    //                         pos = wit->first.find('?', pos + 1);
    //                     }
    //                 }
    // ```
    // A format without rules, or a prefix with a byte that no rule can match (rules are
    // made of `A`-`Z`, `_` and `?`; the prefix has no digit), gives `eAcc_unknown`. LOSAT
    // stops there: the rules and specials themselves are not ported.
    let prefix_can_match = main[..digit_pos]
        .iter()
        .all(|&byte| byte.is_ascii_uppercase() || byte == b'_' || byte == b'?');
    (prefix_can_match && GUIDE_FORMATS.contains(&(digit_pos, digit_count)))
        .then_some(SeqIdLine::GuideFormat)
}

/// What `CSeq_id` makes of the trimmed first line of a record.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2457-2503
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
/// ```
/// A tagged line never reaches the only `Malformatted ID` throw (Seq_id.cpp:2501-2502;
/// `x_Init` throws other texts, such as `Malformatted PDB ID`). A valid local ID is
/// parsed again without `fParse_ValidLocal` (blast_fasta_input.cpp:137-140), which ends
/// at that throw, so a line is FASTA exactly when it is neither tagged, an identified
/// accession, nor a whitelisted general ID.
pub(super) fn classify(line: &[u8]) -> SeqIdLine {
    let line = super::stream::trim_input_space(line);
    // blast_fasta_input.cpp:132: `!line.empty() && isalnum(line.data()[0]&0xff)`.
    if line.first().is_none_or(|&first| !alnum(first)) {
        return SeqIdLine::Fasta;
    }
    if has_fasta_tag(line) {
        return SeqIdLine::FastaTag;
    }
    if let Some(kind) = identify_accession(line) {
        return kind;
    }
    if let Some(colon) = line.iter().position(|&byte| byte == b':') {
        let db = line[..colon].to_ascii_uppercase();
        if RAW_DBTAGS.iter().any(|&known| known == db.as_slice()) {
            return SeqIdLine::General;
        }
    }
    SeqIdLine::Fasta
}

/// LOSAT's rejection of a first line that NCBI may read as a Seq-id, or `None` when NCBI
/// reads it as FASTA (`CSeq_id` throws "Malformatted ID").
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:132-151
/// ```c++
///         const string line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
///         if ( !line.empty() && isalnum(line.data()[0]&0xff) ) {
///             try {
///                 CRef<CSeq_id> id(new CSeq_id(line, (CSeq_id::fParse_AnyRaw |
/// 							CSeq_id::fParse_ValidLocal)));
/// 		if (id->IsLocal()  &&  !NStr::StartsWith(line, "lcl|") ) {
///                     // Expected to throw an exception.
///                     id.Reset(new CSeq_id(line));
/// 		}
///                 CRef<CBioseq> bioseq(x_CreateBioseq(id));
///                 CRef<CSeq_entry> retval(new CSeq_entry());
///                 retval->SetSeq(*bioseq);
///                 return retval;
///             } catch (const CSeqIdException& e) {
///                 if (NStr::Find(e.GetMsg(), "Malformatted ID") != NPOS) {
///                     // This is probably just plain fasta, so just
///                     // defer to CFastaReader
///                 } else {
///                     throw;
///                 }
/// ```
/// Every line that `classify` does not call FASTA is rejected: NCBI fetches it, skips it
/// (a query whose ID is not found, or another rethrown exception) or fails (a subject).
pub(super) fn reject_seq_id_line(line: &[u8], config: &ReaderConfig) -> Option<anyhow::Error> {
    if classify(line) == SeqIdLine::Fasta {
        return None;
    }
    Some(anyhow::anyhow!(
        "the first line of the {role} ({:?}) is not a defline and may be a sequence identifier that NCBI BLAST+ fetches through a data loader (from GenBank or a BLAST database), which is not supported by LOSAT's {program} (start the {role} with a '>' defline, or set DATA_LOADERS=none in the [BLAST] section of .ncbirc)",
        String::from_utf8_lossy(line),
        role = config.role,
        program = config.program,
    ))
}
