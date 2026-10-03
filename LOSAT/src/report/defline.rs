//! The subject title of the pairwise report (outfmt 0) as NCBI makes it for a nucleotide
//! subject read from FASTA: `CDeflineGenerator::GenerateDefline` over the defline, which
//! is the title of a subject read without `-parse_deflines` (its ID included).
//!
//! NCBI reference: c++/src/objmgr/util/create_defline.cpp:3952-3960
//! ```c
//!     if (! m_Reconstruct) {
//!         // x_SetFlags set m_MainTitle from a suitable descriptor, if any;
//!         // now strip trailing periods, commas, semicolons, and spaces.
//!         size_t pos = m_MainTitle.find_last_not_of (".,;~ ");
//!         if (pos != NPOS) {
//!             m_MainTitle.erase (pos + 1);
//!         }
//! ```
//! NCBI reference: c++/src/objmgr/util/create_defline.cpp:4050-4095
//! ```c
//!     if (! (flags & fLeavePrefixSuffix)) {
//!         // remove TPA or TSA prefix, will rely on other data in record to set
//!         for (size_t i = 0; i < sizeof (s_tpaPrefixList) / sizeof (const char*); i++) {
//!             string str = s_tpaPrefixList [i];
//!             if (NStr::StartsWith (m_MainTitle, str, NStr::eNocase)) {
//!                 m_MainTitle.erase (0, str.length());
//!                 // strip leading spaces remaining after removal of old MAG before TPA or TSA prefixes
//!                 m_MainTitle.erase (0, m_MainTitle.find_first_not_of (' '));
//!             }
//!         }
//!     }
//!
//!     // strip leading spaces remaining after removal of old TPA or TSA prefixes
//!     m_MainTitle.erase (0, m_MainTitle.find_first_not_of (' '));
//!
//!     CStringUTF8 decoded = NStr::HtmlDecode (m_MainTitle);
//!
//!     // strip trailing commas, semicolons, and spaces (period may be an sp.
//!     // species)
//!     size_t pos = decoded.find_last_not_of (",;~ ");
//!     if (pos != NPOS) {
//!         decoded.erase (pos + 1);
//!     }
//!     ...
//!     // produce final result
//!     string penult = mag + prefix + decoded + suffix;
//!
//!     x_CleanAndCompress (final, penult, m_IsAA);
//! ```
//! A local FASTA subject has no MolInfo, source or metagenome data, so `mag`, `prefix` and
//! `suffix` are empty; its title is not empty (LOSAT rejects empty deflines), so it is not
//! capitalized. `NStr::HtmlDecode` changes only character references; the callers reject
//! a subject whose title NCBI decodes (`ncbi_nucleotide_title_is_decoded`). The alignment
//! heading removes the prefixes; the description table keeps them (`fLeavePrefixSuffix`,
//! showdefline.cpp:498). For some titles of commas, semicolons, tildes and spaces, NCBI's
//! `x_CleanAndCompress` lets its count of the remaining letters wrap and reads past the
//! end of the string (and crashes); `clean_and_compress` stops at the end of the string
//! there (approved exception 2 of PD-LOSAT-NCBI-DEFECTS).

/// NCBI reference: c++/src/objmgr/util/create_defline.cpp:3431-3446
/// ```c
/// static const char* s_tpaPrefixList [] = {
///   "MAG ",
///   "MAG:",
///   "MULTISPECIES:",
///   "TLS:",
///   "TPA:",
///   "TPA_exp:",
///   "TPA_inf:",
///   "TPA_reasm:",
///   "TPA_asm:",
///   "TSA:",
///   "UNVERIFIED_ORG:",
///   "UNVERIFIED_ASMBLY:",
///   "UNVERIFIED_CONTAM:",
///   "UNVERIFIED:"
/// };
/// ```
const TPA_PREFIXES: [&str; 14] = [
    "MAG ",
    "MAG:",
    "MULTISPECIES:",
    "TLS:",
    "TPA:",
    "TPA_exp:",
    "TPA_inf:",
    "TPA_reasm:",
    "TPA_asm:",
    "TSA:",
    "UNVERIFIED_ORG:",
    "UNVERIFIED_ASMBLY:",
    "UNVERIFIED_CONTAM:",
    "UNVERIFIED:",
];

/// The title of a nucleotide subject with the defline `defline` (the text after `>`),
/// for the alignment heading (`leave_prefix` false) or the description table (true).
pub fn ncbi_nucleotide_title(defline: &str, leave_prefix: bool) -> String {
    let (cleaned, _) = clean_and_compress(
        &nucleotide_title_before_cleanup(defline, leave_prefix),
        false,
    );
    String::from_utf8(cleaned).expect("ASCII deflines")
}

/// Whether NCBI's `x_CleanAndCompress` reads past the end of the outfmt 0 title of a
/// nucleotide subject with the defline `defline`, in the alignment heading or in the
/// description table (NCBI crashes when it writes the title of such a subject with hits).
/// Approved exception 2 of PD-LOSAT-NCBI-DEFECTS covers BLASTN only; the other programs
/// reject such subjects.
pub fn ncbi_nucleotide_title_reads_past_end(defline: &str) -> bool {
    [false, true].into_iter().any(|leave_prefix| {
        clean_and_compress(
            &nucleotide_title_before_cleanup(defline, leave_prefix),
            false,
        )
        .1
    })
}

/// The title of `ncbi_nucleotide_title` before `x_CleanAndCompress`.
fn nucleotide_title_before_cleanup(defline: &str, leave_prefix: bool) -> Vec<u8> {
    let mut title = nucleotide_title_before_decoding(defline, leave_prefix);
    trim_end_of(&mut title, b",;~ ");
    title
}

/// The title that `GenerateDefline` gives `NStr::HtmlDecode` (`m_MainTitle`, see the
/// module), for the alignment heading (`leave_prefix` false) or the description table.
fn nucleotide_title_before_decoding(defline: &str, leave_prefix: bool) -> Vec<u8> {
    let mut title = defline.as_bytes().to_vec();
    trim_end_of(&mut title, b".,;~ ");
    if !leave_prefix {
        for prefix in TPA_PREFIXES {
            if title.len() >= prefix.len()
                && title[..prefix.len()].eq_ignore_ascii_case(prefix.as_bytes())
            {
                title.drain(..prefix.len());
                trim_leading_spaces(&mut title);
            }
        }
    }
    trim_leading_spaces(&mut title);
    title
}

/// Whether `NStr::HtmlDecode` changes the outfmt 0 title of a nucleotide subject with the
/// defline `defline`, in the alignment heading or in the description table. LOSAT does
/// not decode, so the callers reject such a subject where NCBI writes its title.
pub fn ncbi_nucleotide_title_is_decoded(defline: &str) -> bool {
    [false, true].into_iter().any(|leave_prefix| {
        html_decode_changes(&nucleotide_title_before_decoding(defline, leave_prefix))
    })
}

/// Rejects the outfmt 0 title of the subject record `record` (from 1) of `program` that
/// LOSAT does not write as NCBI: a title that NCBI decodes (`NStr::HtmlDecode`), or one
/// that NCBI's `x_CleanAndCompress` reads past (NCBI crashes; approved exception 2 of
/// PD-LOSAT-NCBI-DEFECTS covers BLASTN only). NCBI makes the titles only of the subjects
/// that a report shows, those with hits, so the callers check those subjects after the
/// search.
///
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1540
/// ```c
///         x_DisplayDeflines(aln_set, itr_num, prev_seqids);
/// ```
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1550-1551
/// ```c
///     CSeq_align_set copy_aln_set;
///     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
/// ```
pub fn check_shown_subject_title(
    defline: &str,
    record: usize,
    program: &str,
) -> anyhow::Result<()> {
    if ncbi_nucleotide_title_is_decoded(defline) {
        anyhow::bail!(
            "subject record {record} has an HTML character reference (such as &amp;) in its defline, which NCBI BLAST+ decodes in the outfmt 0 titles; this is not supported by LOSAT's {program}"
        );
    }
    if ncbi_nucleotide_title_reads_past_end(defline) {
        anyhow::bail!(
            "subject record {record} has a defline of punctuation that NCBI BLAST+ reads past its end when it writes the subject's outfmt 0 title (it crashes); this is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}

/// The entity names of `NStr::HtmlDecode`, in NCBI's order.
///
/// NCBI reference: c++/src/corelib/ncbistr.cpp:4223-4230
/// ```c
/// static struct tag_HtmlEntities
/// {
///     TUnicodeSymbol u;
///     const char*    s;
/// }
/// const s_HtmlEntities[] = {
///     {    9, "Tab" },
///     {   10, "NewLine" },
/// ```
/// (and the rows after them, to `{ 9830, "diams" }` at ncbistr.cpp:4508.)
#[rustfmt::skip]
const NCBI_HTML_ENTITY_NAMES: [&str; 280] = [
    "Tab", "NewLine", "excl", "quot", "num", "dollar", "percnt", "amp", "apos", "lpar",
    "rpar", "ast", "plus", "comma", "period", "sol", "colon", "semi", "lt", "equals", "gt",
    "quest", "commat", "lsqb", "bsol", "rsqb", "Hat", "lowbar", "grave", "lcub", "verbar",
    "rcub", "nbsp", "iexcl", "cent", "pound", "curren", "yen", "brvbar", "sect", "uml",
    "copy", "ordf", "laquo", "not", "shy", "reg", "macr", "deg", "plusmn", "sup2", "sup3",
    "acute", "micro", "para", "middot", "cedil", "sup1", "ordm", "raquo", "frac14",
    "frac12", "frac34", "iquest", "Agrave", "Aacute", "Acirc", "Atilde", "Auml", "Aring",
    "AElig", "Ccedil", "Egrave", "Eacute", "Ecirc", "Euml", "Igrave", "Iacute", "Icirc",
    "Iuml", "ETH", "Ntilde", "Ograve", "Oacute", "Ocirc", "Otilde", "Ouml", "times",
    "Oslash", "Ugrave", "Uacute", "Ucirc", "Uuml", "Yacute", "THORN", "szlig", "agrave",
    "aacute", "acirc", "atilde", "auml", "aring", "aelig", "ccedil", "egrave", "eacute",
    "ecirc", "euml", "igrave", "iacute", "icirc", "iuml", "eth", "ntilde", "ograve",
    "oacute", "ocirc", "otilde", "ouml", "divide", "oslash", "ugrave", "uacute", "ucirc",
    "uuml", "yacute", "thorn", "yuml", "OElig", "oelig", "Scaron", "scaron", "Yuml", "fnof",
    "circ", "tilde", "Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", "Theta",
    "Iota", "Kappa", "Lambda", "Mu", "Nu", "Xi", "Omicron", "Pi", "Rho", "Sigma", "Tau",
    "Upsilon", "Phi", "Chi", "Psi", "Omega", "alpha", "beta", "gamma", "delta", "epsilon",
    "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi",
    "rho", "sigmaf", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega", "thetasym",
    "upsih", "piv", "ensp", "emsp", "thinsp", "zwnj", "zwj", "lrm", "rlm", "ndash", "mdash",
    "lsquo", "rsquo", "sbquo", "ldquo", "rdquo", "bdquo", "dagger", "Dagger", "bull",
    "hellip", "permil", "prime", "Prime", "lsaquo", "rsaquo", "oline", "frasl", "euro",
    "weierp", "image", "real", "trade", "alefsym", "larr", "uarr", "rarr", "darr", "harr",
    "crarr", "lArr", "uArr", "rArr", "dArr", "hArr", "forall", "part", "exist", "empty",
    "nabla", "isin", "notin", "ni", "prod", "sum", "minus", "lowast", "radic", "prop",
    "infin", "ang", "and", "or", "cap", "cup", "int", "there4", "sim", "cong", "asymp",
    "ne", "equiv", "le", "ge", "sub", "sup", "nsub", "sube", "supe", "oplus", "otimes",
    "perp", "sdot", "lceil", "rceil", "lfloor", "rfloor", "lang", "rang", "loz", "spades",
    "clubs", "hearts", "diams",
];

/// Whether `NStr::HtmlDecode` decodes a character reference of the ASCII text `text` (its
/// result then differs from `text`), scanning as NCBI scans.
///
/// NCBI reference: c++/src/corelib/ncbistr.cpp:4543-4585
/// ```c
///         if (i != e && ch == '&') {
///             CTempString::const_iterator start_of_entity, end_of_entity, itmp;
///             end_of_entity = itmp = i;
///             bool ent, dec, hex, parsed=false;
///             ent = isalpha((unsigned char)(*itmp)) != 0;
///             dec = !ent && *itmp == '#' && ++itmp != e &&
///                   isdigit((unsigned char)(*itmp)) != 0;
///             hex = !dec && itmp != e &&
///                   (*itmp == 'x' || *itmp == 'X') && ++itmp != e &&
///                    isxdigit((unsigned char)(*itmp)) != 0;
///             start_of_entity = itmp;
///
///             if (itmp != e && (ent || dec || hex)) {
///                 // do not look too far
///                 for (int len=0; len<16 && itmp != e; ++len, ++itmp) {
///                     if (*itmp == '&' || *itmp == '#') {
///                         break;
///                     }
///                     if (*itmp == ';') {
///                         end_of_entity = itmp;
///                         break;
///                     }
///                     ent = ent && isalnum( (unsigned char)(*itmp)) != 0;
///                     dec = dec && isdigit( (unsigned char)(*itmp)) != 0;
///                     hex = hex && isxdigit((unsigned char)(*itmp)) != 0;
///                 }
///                 if (end_of_entity != i && (ent || dec || hex)) {
///                     uch = 0;
///                     if (ent) {
///                         string entity(start_of_entity, end_of_entity);
///                         const struct tag_HtmlEntities* p = s_HtmlEntities;
///                         for ( ; p->u != 0; ++p) {
///                             if (entity.compare(p->s) == 0) {
///                                 uch = p->u;
///                                 parsed = true;
///                                 result |= fHtmlDec_CharRef_Entity;
///                                 break;
///                             }
///                         }
///                     } else {
///                         parsed = true;
/// ```
/// A numeric reference is always decoded. As in NCBI, an `x` before a hexadecimal digit
/// moves the start of a named entity past the `x`, so `&xi;` is not decoded.
fn html_decode_changes(text: &[u8]) -> bool {
    let e = text.len();
    let mut i = 0;
    while i < e {
        let ch = text[i];
        i += 1;
        if i == e || ch != b'&' {
            continue;
        }
        let mut itmp = i;
        let mut end_of_entity = i;
        let mut ent = text[itmp].is_ascii_alphabetic();
        let mut dec = false;
        if !ent && text[itmp] == b'#' {
            itmp += 1;
            dec = itmp != e && text[itmp].is_ascii_digit();
        }
        let mut hex = false;
        if !dec && itmp != e && matches!(text[itmp], b'x' | b'X') {
            itmp += 1;
            hex = itmp != e && text[itmp].is_ascii_hexdigit();
        }
        let start_of_entity = itmp;
        if itmp == e || !(ent || dec || hex) {
            continue;
        }
        let mut len = 0;
        while len < 16 && itmp != e {
            let c = text[itmp];
            if c == b'&' || c == b'#' {
                break;
            }
            if c == b';' {
                end_of_entity = itmp;
                break;
            }
            ent = ent && c.is_ascii_alphanumeric();
            dec = dec && c.is_ascii_digit();
            hex = hex && c.is_ascii_hexdigit();
            len += 1;
            itmp += 1;
        }
        if end_of_entity != i && (ent || dec || hex) {
            if !ent {
                return true;
            }
            let entity = &text[start_of_entity..end_of_entity];
            if NCBI_HTML_ENTITY_NAMES
                .iter()
                .any(|name| name.as_bytes() == entity)
            {
                return true;
            }
        }
    }
    false
}

/// `find_last_not_of(chars)` then `erase(pos + 1)`; a text made only of `chars` stays.
fn trim_end_of(text: &mut Vec<u8>, chars: &[u8]) {
    if let Some(pos) = text.iter().rposition(|byte| !chars.contains(byte)) {
        text.truncate(pos + 1);
    }
}

/// `erase(0, find_first_not_of(' '))`: a text of spaces only becomes empty.
fn trim_leading_spaces(text: &mut Vec<u8>) {
    let first = text
        .iter()
        .position(|&byte| byte != b' ')
        .unwrap_or(text.len());
    text.drain(..first);
}

/// NCBI's `x_CleanAndCompress`, byte for byte (the titles are ASCII), except that it stops
/// at the end of the string where NCBI's count of the remaining letters wraps (see below);
/// the flag is whether it wrapped.
///
/// NCBI reference: c++/src/objmgr/util/create_defline.cpp:219-312
/// ```c
///     char curr = *in++; // initialize with first character
///     left--;
///
///     char next = 0;
///     Uint2 two_chars = curr; // this is two bytes storage where we see current and previous symbols
///
///     while (left > 0) {
///         next = *in++;
///
///         two_chars = Uint2((two_chars << 8) | next);
///
///         switch (two_chars)
///         {
///         case twocommas: // replace double commas with comma+space
///             *out++ = curr;
///             next = ' ';
///             break;
///         case twospaces: // skip multispaces (only print last one)
///             break;
///         case bracket_space: // skip space after bracket
///             next = curr;
///             two_chars = curr;
///             break;
///         case space_bracket: // skip space before bracket
///             break;
///         case space_comma:
///         case space_semicolon: // swap characters
///             *out++ = next;
///             next = curr;
///             two_chars = curr;
///             break;
///         case comma_space:
///             *out++ = curr;
///             *out++ = ' ';
///             while (next == ' ' || next == ',') {
///                 next = *in;
///                 in++;
///                 left--;
///             }
///             two_chars = next;
///             break;
///         case semicolon_space:
///     ...
///         default:
///             *out++ = curr;
///             break;
///         }
///
///         curr = next;
///         left--;
///     }
///
///     if (curr > 0 && curr != ' ') {
///         *out++ = curr;
///     }
/// ```
fn clean_and_compress(input: &[u8], is_protein: bool) -> (Vec<u8>, bool) {
    let mut start = 0;
    let mut end = input.len();
    while start < end && input[start] == b' ' {
        start += 1;
    }
    while end > start && input[end - 1] == b' ' {
        end -= 1;
    }
    let text = &input[start..end];
    let mut out = Vec::with_capacity(text.len());
    if text.is_empty() {
        return (out, false);
    }
    // Past the end reads the terminating NUL of the C++ string.
    let at = |index: usize| text.get(index).copied().unwrap_or(0);
    let mut index = 0;
    let mut left = text.len();
    let mut wrapped = false;
    let mut curr = at(index);
    index += 1;
    left -= 1;
    let mut two_chars = u16::from(curr);
    while left > 0 {
        let mut next = at(index);
        index += 1;
        two_chars = (two_chars << 8) | u16::from(next);
        match two_chars.to_be_bytes() {
            [b',', b','] => {
                out.push(curr);
                next = b' ';
            }
            [b' ', b' '] | [b' ', b')'] => {}
            [b'(', b' '] => {
                next = curr;
                two_chars = u16::from(curr);
            }
            [b' ', b','] | [b' ', b';'] => {
                out.push(next);
                next = curr;
                two_chars = u16::from(curr);
            }
            [separator @ (b',' | b';'), b' '] => {
                out.push(curr);
                out.push(b' ');
                while next == b' ' || next == separator {
                    next = at(index);
                    index += 1;
                    // Where NCBI's `left` (a size_t) wraps, for a run of spaces and
                    // separators that reaches the end of the string, its loop runs past the
                    // string (it reads the terminating NUL and on, and crashes). LOSAT stops
                    // at the end of the string (approved exception 2 of
                    // PD-LOSAT-NCBI-DEFECTS): `next` is the NUL and nothing more is written.
                    wrapped |= left == 0;
                    left = left.saturating_sub(1);
                }
                two_chars = u16::from(next);
            }
            _ => out.push(curr),
        }
        curr = next;
        wrapped |= left == 0;
        left = left.saturating_sub(1);
    }
    if curr > 0 && curr != b' ' {
        out.push(curr);
    }
    if is_protein {
        let replaced = String::from_utf8_lossy(&out)
            .replace(". [", " [")
            .replace(", [", " [");
        return (replaced.into_bytes(), wrapped);
    }
    (out, wrapped)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn titles_follow_ncbi_defline_generator() {
        // NCBI BLAST+ 2.17.0 outfmt 0 (alignment heading, description table).
        for (defline, heading, description) in [
            ("s1 abc.", "s1 abc", "s1 abc"),
            ("s1 abc ,def", "s1 abc, def", "s1 abc, def"),
            ("s1 a  b", "s1 a b", "s1 a b"),
            ("s1 (a )", "s1 (a)", "s1 (a)"),
            ("s1,,x", "s1, x", "s1, x"),
            ("TPA: s1 x", "s1 x", "TPA: s1 x"),
            ("s1 TPA: x", "s1 TPA: x", "s1 TPA: x"),
            ("s1 MAG x", "s1 MAG x", "s1 MAG x"),
            ("s1 x ; y", "s1 x; y", "s1 x; y"),
            ("s1 x,  y", "s1 x, y", "s1 x, y"),
            ("s1 x~", "s1 x", "s1 x"),
            ("s1 E. coli sp.", "s1 E. coli sp", "s1 E. coli sp"),
            ("s1 a, ,b", "s1 a, b", "s1 a, b"),
            ("s1 x;;y", "s1 x;;y", "s1 x;;y"),
            ("s1  lead", "s1 lead", "s1 lead"),
            ("s1 a ( b", "s1 a (b", "s1 a (b"),
            ("s1..x", "s1..x", "s1..x"),
            ("s1.x.", "s1.x", "s1.x"),
        ] {
            assert_eq!(
                ncbi_nucleotide_title(defline, false),
                heading,
                "{defline:?}"
            );
            assert_eq!(
                ncbi_nucleotide_title(defline, true),
                description,
                "{defline:?}"
            );
        }
        // NCBI BLAST+ 2.17.0 crashes on the outfmt 0 titles of the first group
        // (x_CleanAndCompress reads past the string); LOSAT stops at the end of the string
        // (approved exception 2 of PD-LOSAT-NCBI-DEFECTS). The second group is NCBI's.
        for (defline, title) in [
            (", ,", ", "),
            ("; ;", "; "),
            ("~, ,", "~, "),
            (",, ,", ",  "),
            (", ,,", ", "),
            (";  ;", "; "),
            (", , ,", ", "),
            (",~, ,", ",~, "),
            (", ;", ", ;"),
            (",,", ","),
            ("a, ,b", "a, b"),
            (", ,a", ", a"),
        ] {
            assert_eq!(ncbi_nucleotide_title(defline, false), title, "{defline:?}");
        }
        for defline in [
            ", ,", "; ;", "~, ,", ",, ,", ", ,,", ";  ;", ", , ,", ",~, ,",
        ] {
            assert!(ncbi_nucleotide_title_reads_past_end(defline), "{defline:?}");
        }
        for defline in [", ;", ",,", "a, ,b", ", ,a", "id, ,", "MAG: , ,x"] {
            assert!(
                !ncbi_nucleotide_title_reads_past_end(defline),
                "{defline:?}"
            );
        }
    }

    // NCBI reference: c++/src/corelib/ncbistr.cpp:4543-4590 (NStr::HtmlDecode): table
    // entities and numeric references only, in the title after its final `.,;~ ` are
    // trimmed (create_defline.cpp:3952-3960), with NCBI's start of an entity after `x`.
    #[test]
    fn titles_that_ncbi_decodes_are_found() {
        for defline in [
            "a&amp;b",
            "a&#38;b",
            "a&#x26;b",
            "a&#X26;b",
            "x &lt; y",
            "s &Tab;t",
            "a&#0;b",
            "TPA: &amp;x",
            "a &amp;; b",
        ] {
            assert!(ncbi_nucleotide_title_is_decoded(defline), "{defline}");
        }
        for defline in [
            "a & b",
            "a&b",
            "a&;",
            "a&#;",
            "a&amp b",
            "a&x#y;",
            "R&D; x",
            "a&foo;b",
            "a&amp;",
            "a&amp;;",
            "s &xi;t",
            "a&X41;b",
            "a&amp#;b",
            "a&ampxxxxxxxxxxxxxxxx;b",
            "a&#12345678901234567;b",
        ] {
            assert!(!ncbi_nucleotide_title_is_decoded(defline), "{defline}");
        }
    }
}
