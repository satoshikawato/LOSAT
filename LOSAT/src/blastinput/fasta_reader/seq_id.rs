//! The first line of a record that `CBlastInputReader` tries as a Seq-id.
//!
//! NCBI reads such a line as the ID of a sequence to fetch with its data loaders (GenBank
//! over the network, or a BLAST database), which LOSAT does not do; it is an explicit
//! rejection. A line that NCBI's `CSeq_id` parser rejects as a malformatted ID is read as
//! FASTA data, as NCBI does.
//!
//! The parse is `CSeq_id::Set` (Seq_id.cpp:2457-2551) with `fParse_AnyRaw`
//! and `fParse_ValidLocal`, then the retry with `fParse_AnyRaw` when the result is a local ID.
//! Every branch is ported except the accession guide (`SAccGuide::Find` over the 1458
//! prefix rules and the special ranges of `accguide2.inc`): a line in one of the guide's
//! formats is rejected without the lookup (`SeqIdLine::GuideFormat`), which also rejects
//! some lines that the guide does not know and NCBI reads as FASTA (`ZZ123456`).

/// What the `CSeq_id` constructor with `fParse_AnyRaw` and `fParse_ValidLocal`, and its retry,
/// make of a line.
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
// (commented-out aliases omitted). The tags that give a type other than `e_not_set` and
// that a 2- or 3-byte prefix can match; `gibbsq` is six bytes.
const FASTA_TAGS: [&[u8]; 22] = [
    b"bbm", b"bbs", b"dbj", b"emb", b"gb", b"gi", b"gim", b"gnl", b"gpp", b"lcl", b"nat", b"pat",
    b"pdb", b"pgp", b"pir", b"prf", b"ref", b"sp", b"tpd", b"tpe", b"tpg", b"tr",
];

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:96-124
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
/// The flags are `fParse_AnyRaw`, `fParse_ValidLocal` and `fParse_FallbackOK`
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
            #[expect(
                clippy::if_same_then_else,
                reason = "NCBI returns eAcc_pdb both for a delimited ID and in its switch's case 4 (Seq_id.cpp:1683-1710)"
            )]
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
        let oqp = matches!(main[0], b'O' | b'P' | b'Q');
        #[expect(
            clippy::if_same_then_else,
            reason = "NCBI has the two UniProt patterns as two branches that both return eAcc_swissprot (Seq_id.cpp:1713-1727)"
        )]
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
    if digit_pos == 0 {
        return (!has_version && main[0] != b'0' && main.iter().all(|&byte| digit(byte)))
            .then_some(SeqIdLine::Gi);
    }
    let digit_count = main_size - digit_pos - usize::from(flag_char);
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1370-1398
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
/// A tagged line never reaches the only `Malformatted ID` throw (Seq_id.cpp:2501-2502;
/// `x_Init` throws other texts, such as `Malformatted PDB ID`). A valid local ID is
/// parsed again without `fParse_ValidLocal` (blast_fasta_input.cpp:137-140), which ends
/// at that throw, so a line is FASTA exactly when it is neither tagged, an identified
/// accession, nor a whitelisted general ID.
pub(super) fn classify(line: &[u8]) -> SeqIdLine {
    // `NStr::TruncateSpaces_Unsafe` of `the_id_in` with `NStr::eTrunc_Both` (Seq_id.cpp:2459).
    let line = crate::ncbi_str::truncate_spaces(line);
    // blast_fasta_input.cpp:132: the line is not empty and its first byte (masked to 8 bits) is
    // alphanumeric (`isalnum`).
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
        if RAW_DBTAGS.contains(&db.as_slice()) {
            return SeqIdLine::General;
        }
    }
    SeqIdLine::Fasta
}

/// LOSAT's rejection of a first line that NCBI may read as a Seq-id, or `None` when NCBI
/// reads it as FASTA (`CSeq_id` throws "Malformatted ID").
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:132-151
/// Every line that `classify` does not call FASTA is rejected: NCBI fetches it, skips it
/// (a query whose ID is not found, or another rethrown exception) or fails (a subject).
/// `line` is the line as read (`classify` trims it); `program` and `role` name the reader in
/// the rejection. Both readers ask here: the shared reader (`reader::read_one_seq`, with its
/// data loaders on), BLASTX's too (`algorithm::blastx::input_source`). BLASTX's web ABI v1
/// (`algorithm::blastx::web`, frozen) reads without the data loaders: it reads such a line as data,
/// as before.
pub(crate) fn reject_seq_id_line(line: &[u8], program: &str, role: &str) -> Option<anyhow::Error> {
    if classify(line) == SeqIdLine::Fasta {
        return None;
    }
    Some(anyhow::anyhow!(
        "the first line of the {role} ({:?}) is not a defline and may be a sequence identifier that NCBI BLAST+ fetches through a data loader (from GenBank or a BLAST database), which is not supported by LOSAT's {program} (start the {role} with a '>' defline, or set DATA_LOADERS=none in the [BLAST] section of .ncbirc)",
        String::from_utf8_lossy(crate::ncbi_str::truncate_spaces(line)),
    ))
}
