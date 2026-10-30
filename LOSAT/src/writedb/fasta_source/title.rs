//! The title that `CWriteDB::ExtractBioseqDeflines` gives the Blast-def-line of a FASTA
//! record without `-parse_seqids`: `CFastaReader` keeps the definition line in a user
//! object as `NStr::PrintableString` writes it, and `x_GetFastaReaderDeflines` reads it
//! back (`c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp`, `c++/src/corelib/ncbistr.cpp`).

use std::borrow::Cow;

/// `NStr::PrintableString(str)` with its default mode (`fNewLine_Quote |
/// fNonAscii_Passthru`, which are 0) and `eLanguage_C`.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3732-3845
/// ```c++
/// static string s_PrintableString(const CTempString    str,
///                                 NStr::TPrintableMode mode,
///                                 ELanguage            lang)
/// {
///     unique_ptr<CNcbiOstrstream> out;
///     SIZE_TYPE i, j = 0;
///
///     for (i = 0;  i < str.size();  ++i) {
///         bool octal = false;
///         char c = str[i];
///         switch (c) {
///         case '\a':
///             if (lang == eLanguage_C)
///                 c = 'a';
///             else
///                 octal = true;
///             break;
///         case '\b':
///             c = 'b';
///             break;
///         case '\f':
///             c = 'f';
///             break;
///         case '\r':
///             c = 'r';
///             break;
///         case '\t':
///             c = 't';
///             break;
///         case '\v':
///             c = 'v';
///             break;
///         case '\n':
///             if (!(mode & NStr::fNewLine_Passthru))
///                 c = 'n';
///             /*FALLTHRU*/
///         case '\\':
///         case '\'':
///         case '"':
///             break;
///         case '&':
///             if (lang == eLanguage_Javascript)
///                 break;
///             continue;
///         case '?':
///             if (lang == eLanguage_C) {
///                 if (i  &&  str[i - 1] == '?')
///                     break;
///                 if (i < str.size() - 1  &&  str[i + 1] == '?')
///                     break;
///             }
///             continue;
///         default:
///             if (!isascii((unsigned char) c)) {
///                 if (mode & NStr::fNonAscii_Quote) {
///                     octal = true;
///                     break;
///                 }
///             }
///             if (!isprint((unsigned char) c)) {
///                 octal = true;
///                 break;
///             }
///             continue;
///         }
///         if (!out.get()) {
///             out.reset(new CNcbiOstrstream);
///         }
///         if (i > j) {
///             out->write(str.data() + j, i - j);
///         }
///         out->put('\\');
///         if (c == '\n') {
///             out->write("n\\\n", 3);
///         } else if (octal) {
///             bool reduce;
///             if (!(mode & NStr::fPrintable_Full)) {
///                 reduce = (i == str.size() - 1  ||
///                           str[i + 1] < '0' || '7' < str[i + 1] ? true : false);
///             } else {
///                 reduce = false;
///             }
///             unsigned char v;
///             char val[3];
///             int k = 0;
///             v =  (unsigned char) c >> 6;
///             if (v  ||  !reduce) {
///                 val[k++] = char('0' + v);
///                 reduce = false;
///             }
///             v = ((unsigned char) c >> 3) & 7;
///             if (v  ||  !reduce) {
///                 val[k++] = char('0' + v);
///             }
///             v =  (unsigned char) c       & 7;
///             val[k++] = char('0' + v);
///             out->write(val, k);
///         } else {
///             out->put(c);
///         }
///         j = i + 1;
///     }
///     if (j  &&  i > j) {
///         _ASSERT(out.get());
///         out->write(str.data() + j, i - j);
///     }
///     if (out.get()) {
///         // Return encoded string
///         return CNcbiOstrstreamToString(*out);
///     }
///
///     // All characters are good - return (a copy of) the original string
///     return str;
/// }
/// ```
/// `isprint` is the "C" locale's (0x20 to 0x7E); a locale that prints more bytes would
/// leave them as they are, and `parse_escapes` reads either form as the same byte. The
/// `c == '\n'` branch after the backslash needs `fNewLine_Passthru` (without it `c` is
/// already `n`).
pub(super) fn printable_string(str: &[u8]) -> Cow<'_, [u8]> {
    let mut out: Option<Vec<u8>> = None;
    let mut j = 0;
    for (i, &byte) in str.iter().enumerate() {
        let mut octal = false;
        let mut c = byte;
        match byte {
            0x07 => c = b'a',
            0x08 => c = b'b',
            0x0c => c = b'f',
            b'\r' => c = b'r',
            b'\t' => c = b't',
            0x0b => c = b'v',
            b'\n' => c = b'n',
            b'\\' | b'\'' | b'"' => {}
            b'&' => continue,
            b'?' => {
                let after_one = i > 0 && str[i - 1] == b'?';
                let before_one = i + 1 < str.len() && str[i + 1] == b'?';
                if !after_one && !before_one {
                    continue;
                }
            }
            0x20..=0x7e => continue,
            _ => octal = true,
        }
        let out = out.get_or_insert_with(|| Vec::with_capacity(str.len() + 8));
        out.extend_from_slice(&str[j..i]);
        out.push(b'\\');
        if octal {
            let mut reduce = i == str.len() - 1 || !(b'0'..=b'7').contains(&str[i + 1]);
            let v = c >> 6;
            if v != 0 || !reduce {
                out.push(b'0' + v);
                reduce = false;
            }
            let v = (c >> 3) & 7;
            if v != 0 || !reduce {
                out.push(b'0' + v);
            }
            out.push(b'0' + (c & 7));
        } else {
            out.push(c);
        }
        j = i + 1;
    }
    match out {
        Some(mut out) => {
            out.extend_from_slice(&str[j..]);
            Cow::Owned(out)
        }
        None => Cow::Borrowed(str),
    }
}

/// `NStr::Replace` and `NStr::ReplaceInPlace`: every occurrence of `search`, from the
/// left, without reading the replacement again.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3371-3385
/// ```c++
///     } else {
///         // Replacing string is shorter or have the same length.
///         // ReplaceInPlace() can be faster on some platform, but not much,
///         // so we use regular algorithm even for equal lengths here.
///         dst = src;
///         for (SIZE_TYPE count = 0; !(max_replace && count >= max_replace); count++){
///             start_pos = dst.find(search, start_pos);
///             if (start_pos == NPOS)
///                 break;
///             dst.replace(start_pos, search.size(), replace);
///             start_pos += replace.size();
///             if (num_replace)
///                 (*num_replace)++;
///         }
///     }
/// ```
fn replace(src: &[u8], search: &[u8], replacement: &[u8]) -> Vec<u8> {
    let mut dst = Vec::with_capacity(src.len());
    let mut pos = 0;
    while let Some(found) = src[pos..]
        .windows(search.len())
        .position(|window| window == search)
    {
        dst.extend_from_slice(&src[pos..pos + found]);
        dst.extend_from_slice(replacement);
        pos += found + search.len();
    }
    dst.extend_from_slice(&src[pos..]);
    dst
}

/// `NStr::ParseEscapes(str)` with its default mode (`eEscSeqRange_Standard`). An error is
/// the message of NCBI's `CStringException`; the text that `fasta_reader_title` gives it
/// has no escape that fails (see there).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4789-4892
/// ```c++
/// string NStr::ParseEscapes(const CTempString str, EEscSeqRange mode, char user_char)
/// {
///     string out;
///     out.reserve(str.size());  // result string can only be smaller
///     SIZE_TYPE pos = 0;
///     bool is_error = false;
///
///     while (pos < str.size()  ||  !is_error) {
///         SIZE_TYPE pos2 = str.find('\\', pos);
///         if (pos2 == NPOS) {
///             //~ out += str.substr(pos);
///             CTempString sub(str, pos);
///             out += sub;
///             break;
///         }
///         //~ out += str.substr(pos, pos2 - pos);
///         CTempString sub(str, pos, pos2-pos);
///         out += sub;
///         if (++pos2 == str.size()) {
///             NCBI_THROW2(CStringException, eFormat,
///                         "Unterminated escape sequence", pos2);
///         }
///         switch (str[pos2]) {
///         case 'a':  out += '\a';  break;
///         case 'b':  out += '\b';  break;
///         case 'f':  out += '\f';  break;
///         case 'n':  out += '\n';  break;
///         case 'r':  out += '\r';  break;
///         case 't':  out += '\t';  break;
///         case 'v':  out += '\v';  break;
///         case 'x':
///             {{
///                 pos = ++pos2;
///                 while (pos < str.size()
///                        &&  isxdigit((unsigned char) str[pos])) {
///                     pos++;
///                 }
///                 if (pos > pos2) {
///                     SIZE_TYPE len = pos-pos2;
///                     if ((mode == eEscSeqRange_FirstByte) && (len > 2)) {
///                         // Take only 2 first hex-digits
///                         len = 2;
///                         pos = pos2 + 2;
///                     }
///                     unsigned int value =
///                         StringToUInt(CTempString(str, pos2, len), 0, 16);
///                     if ((mode != eEscSeqRange_Standard)  &&  (value > 255)) {
///                         // eEscSeqRange_Standard -- by default
///                         switch (mode) {
///                         case eEscSeqRange_FirstByte:
///                             // Already have right value
///                             break;
///                         case eEscSeqRange_Throw:
///                             NCBI_THROW2(CStringException, eFormat,
///                                 "Escape sequence '" + NStr::PrintableString(CTempString(str, pos2, len)) +
///                                 "' is out of range [0-255]", pos2);
///                             break;
///                         case eEscSeqRange_Errno:
///                             CNcbiError::SetErrno(errno = ERANGE, str);
///                             is_error = true;
///                             continue;
///                         case eEscSeqRange_User:
///                             value = (unsigned)user_char;
///                             break;
///                         default:
///                             NCBI_THROW2(CStringException, eFormat, "Wrong set of flags", pos2);
///                         }
///                     }
///                     out += static_cast<char>(value);
///                 } else {
///                     NCBI_THROW2(CStringException, eFormat,
///                                 "\\x followed by no hexadecimal digits", pos);
///                 }
///             }}
///             continue;
///         case '0':  case '1':  case '2':  case '3':
///         case '4':  case '5':  case '6':  case '7':
///             {{
///                 pos = pos2;
///                 unsigned char c = (unsigned char)(str[pos++] - '0');
///                 while (pos < pos2 + 3  &&  pos < str.size()
///                        &&  str[pos] >= '0'  &&  str[pos] <= '7') {
///                     c = (unsigned char)((c << 3) | (str[pos++] - '0'));
///                 }
///                 out += c;
///             }}
///             continue;
///         case '\n':
///             // quoted EOL means no EOL
///             break;
///         default:
///             out += str[pos2];
///             break;
///         }
///         pos = pos2 + 1;
///     }
///     if (mode == eEscSeqRange_Errno) {
///         if (is_error) {
///             return kEmptyStr;
///         }
///         errno = 0;
///     }
///     return out;
/// }
/// ```
fn parse_escapes(str: &[u8]) -> Result<Vec<u8>, &'static str> {
    let mut out = Vec::with_capacity(str.len());
    let mut pos = 0;
    loop {
        let Some(found) = str[pos..].iter().position(|&byte| byte == b'\\') else {
            out.extend_from_slice(&str[pos..]);
            break;
        };
        let mut pos2 = pos + found;
        out.extend_from_slice(&str[pos..pos2]);
        pos2 += 1;
        if pos2 == str.len() {
            return Err("Unterminated escape sequence");
        }
        match str[pos2] {
            b'a' => out.push(0x07),
            b'b' => out.push(0x08),
            b'f' => out.push(0x0c),
            b'n' => out.push(b'\n'),
            b'r' => out.push(b'\r'),
            b't' => out.push(b'\t'),
            b'v' => out.push(0x0b),
            b'x' => {
                pos2 += 1;
                pos = pos2;
                while pos < str.len() && str[pos].is_ascii_hexdigit() {
                    pos += 1;
                }
                if pos == pos2 {
                    return Err("\\x followed by no hexadecimal digits");
                }
                // `StringToUInt(..., 0, 16)` throws for a value above 4,294,967,295.
                let mut value: u32 = 0;
                for &digit in &str[pos2..pos] {
                    let digit = (digit as char).to_digit(16).unwrap_or(0);
                    value = value
                        .checked_mul(16)
                        .and_then(|value| value.checked_add(digit))
                        .ok_or("Cannot convert string to unsigned int")?;
                }
                out.push(value as u8);
                continue;
            }
            b'0'..=b'7' => {
                pos = pos2;
                let mut c = str[pos] - b'0';
                pos += 1;
                while pos < pos2 + 3 && pos < str.len() && (b'0'..=b'7').contains(&str[pos]) {
                    c = (c << 3) | (str[pos] - b'0');
                    pos += 1;
                }
                out.push(c);
                continue;
            }
            b'\n' => {}
            other => out.push(other),
        }
        pos = pos2 + 1;
    }
    Ok(out)
}

/// The `fasta` string of `x_GetFastaReaderDeflines`: the `DefLine` field with `\t` replaced
/// by three spaces, then its escapes read (`fasta_reader_title` quotes the lines).
pub(crate) fn fasta_reader_fasta(defline_field: &[u8]) -> Result<Vec<u8>, &'static str> {
    parse_escapes(&replace(defline_field, b"\\t", TAB_REPLACEMENT))
}

/// `TAB_REPLACEMENT`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1714
/// ```c++
/// #define TAB_REPLACEMENT "   "
/// ```
const TAB_REPLACEMENT: &[u8] = b"   ";

/// The title of `x_GetFastaReaderDeflines` without `parse_ids`, from the `DefLine` field
/// of the `CFastaReader` user object (`printable_string` of the definition line with its
/// `>`); `None` when the field is empty (no definition line from the user object).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1752-1760,1766-1767,1782-1798
/// ```c++
///                     if (f.CanGetLabel() &&
///                         f.GetLabel().IsStr() &&
///                         f.GetLabel().GetStr() == "DefLine" &&
///                         f.CanGetData() &&
///                         f.GetData().IsStr()) {
///                         fasta = NStr::Replace(f.GetData().GetStr(), "\\t", TAB_REPLACEMENT);
///                         fasta = NStr::ParseEscapes(fasta);
///                         break;
///                     }
/// ...
///     if (fasta.empty())
///         return;
/// ...
///     if (!parse_ids) {
///
///         // Generate an BL_ORD_ID in case no parse is needed
///         CRef<CSeq_id> gnl_id(new CSeq_id());
///         gnl_id->SetGeneral().SetDb("BL_ORD_ID");
///         gnl_id->SetGeneral().SetTag().SetId(0);  // will be filled later
///
///         // Build the local defline.
///         defline.Reset(new CBlast_def_line);
///         defline->SetSeqid().push_back(gnl_id);
///
///         string title(fasta, 1, fasta.size());
///         // Replace ^A with space
///         NStr::ReplaceInPlace(title, "\001", " ");
///         // Replace tabs with three spaces
///         NStr::ReplaceInPlace(title, "\t", TAB_REPLACEMENT);
///         defline->SetTitle(title);
/// ```
/// The field is searched for the two characters `\t` before its escapes are read, so a
/// tab of the definition line (written `\t`) becomes three spaces, and so does a
/// backslash followed by `t` (written `\\t`): the search finds the second backslash with
/// the `t`, and the first backslash then escapes the first of the three spaces. A ctrl-A
/// becomes one space. Every other byte comes back as it was read.
///
/// No escape fails: a backslash that starts an escape is one that `printable_string`
/// wrote, followed by a letter, a quote, `?`, a backslash or octal digits, or it is the
/// backslash left of a replaced `\\t`, followed by a space.
pub(super) fn fasta_reader_title(defline_field: &[u8]) -> Result<Option<Vec<u8>>, &'static str> {
    let fasta = fasta_reader_fasta(defline_field)?;
    if fasta.is_empty() {
        return Ok(None);
    }
    let title = replace(&fasta[1..], b"\x01", b" ");
    Ok(Some(replace(&title, b"\t", TAB_REPLACEMENT)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn title(defline: &[u8]) -> Vec<u8> {
        fasta_reader_title(&printable_string(defline))
            .unwrap()
            .unwrap()
    }

    #[test]
    fn printable_string_escapes() {
        assert_eq!(&*printable_string(b">plain title 1"), b">plain title 1");
        assert!(matches!(printable_string(b">a&b?c"), Cow::Borrowed(_)));
        assert_eq!(&*printable_string(b">a\tb"), b">a\\tb");
        assert_eq!(&*printable_string(b">a\\b"), b">a\\\\b");
        assert_eq!(&*printable_string(b">a'b\"c"), b">a\\'b\\\"c");
        assert_eq!(&*printable_string(b">a??b?"), b">a\\?\\?b?");
        assert_eq!(
            &*printable_string(b"\x07\x08\x0c\r\x0b\n"),
            b"\\a\\b\\f\\r\\v\\n"
        );
        // Octal: as few digits as the next byte allows.
        assert_eq!(&*printable_string(b"\x01z"), b"\\1z");
        assert_eq!(&*printable_string(b"\x011"), b"\\0011");
        assert_eq!(&*printable_string(b"\x00"), b"\\0");
        assert_eq!(&*printable_string(b"\x1f8"), b"\\378");
        assert_eq!(&*printable_string(b"\x7f\x80\xff"), b"\\177\\200\\377");
        assert_eq!(&*printable_string(b"\x09\x0a"), b"\\t\\n");
    }

    #[test]
    fn parse_escapes_reads_printable_string_back() {
        let all: Vec<u8> = (0..=255u8).chain([b'1', b'?', b'?', 0, b'7']).collect();
        assert_eq!(parse_escapes(&printable_string(&all)).unwrap(), all);
        assert_eq!(parse_escapes(b"a\\x41\\x4a4bz").unwrap(), b"aA\x4bz");
        assert_eq!(parse_escapes(b"\\101\\1012\\8\\q").unwrap(), b"AA28q");
        assert_eq!(parse_escapes(b"a\\\nb").unwrap(), b"ab");
        assert_eq!(parse_escapes(b"a\\"), Err("Unterminated escape sequence"));
        assert_eq!(
            parse_escapes(b"a\\xg"),
            Err("\\x followed by no hexadecimal digits")
        );
    }

    /// Titles as makeblastdb stores them (the header file of the databases of the probe
    /// cases `p_title_bytes`, `p_title_spaces`, `p_ctrl_a`, `fs_p_title_*`; a byte that
    /// the ASN.1 writer replaces by `#` is the input's byte here).
    #[test]
    fn titles_of_makeblastdb() {
        assert_eq!(title(b">a t"), b"a t");
        assert_eq!(title(b">"), b"");
        assert_eq!(title(b">  a   two  spaces   "), b"  a   two  spaces   ");
        assert_eq!(title(b">b\t\t"), b"b      ");
        assert_eq!(title(b">a first\x01b second"), b"a first b second");
        assert_eq!(title(b">c x >d y"), b"c x >d y");
        assert_eq!(
            title(b">a\ttab \x01ctl \xc3\xa9 utf8 \x80hi \x7f del"),
            b"a   tab  ctl \xc3\xa9 utf8 \x80hi \x7f del"
        );
        assert_eq!(
            title(b">sp|P12345|ABC_HUMAN desc"),
            b"sp|P12345|ABC_HUMAN desc"
        );
        assert_eq!(title(b">>a"), b">a");
        assert_eq!(title(b">\ta\t\tb\t"), b"   a      b   ");
        assert_eq!(title(b">a\x01\tb"), b"a    b");
        assert_eq!(title(b">b\x01\x01\x01"), b"b   ");
        assert_eq!(title(b">a\x011\x027\x008\x1f0z"), b"a 1\x027\x008\x1f0z");
        for byte in (0..=8u8).chain([0x0b, 0x0c]).chain(0x0e..0x20) {
            let expected = if byte == 1 { b' ' } else { byte };
            assert_eq!(title(&[b'>', b'a', byte, b'z']), [b'a', expected, b'z']);
        }
        // A backslash and `t`.
        assert_eq!(title(b">a\\tb"), b"a   b");
        assert_eq!(title(b">a\\\\tb"), b"a\\   b");
        assert_eq!(title(b">a\\\tb"), b"a\\   b");
        assert_eq!(title(b">\\t"), b"   ");
        assert_eq!(title(b">t\\t\\t"), b"t      ");
        // Other backslashes stay.
        assert_eq!(title(b">a\\nb"), b"a\\nb");
        assert_eq!(title(b">end\\"), b"end\\");
        assert_eq!(title(b">a\\x41b"), b"a\\x41b");
        assert_eq!(title(b">a\\101b"), b"a\\101b");
        assert_eq!(title(b">\\\\"), b"\\\\");
        assert_eq!(title(b">a\\b\\a\\f\\r\\v"), b"a\\b\\a\\f\\r\\v");
        assert_eq!(title(b">a'b\"c"), b"a'b\"c");
        assert_eq!(title(b">a???b?"), b"a???b?");
        assert_eq!(title(b">a&b<c>"), b"a&b<c>");
    }

    #[test]
    fn empty_field_gives_no_title() {
        assert_eq!(fasta_reader_title(b""), Ok(None));
    }
}
