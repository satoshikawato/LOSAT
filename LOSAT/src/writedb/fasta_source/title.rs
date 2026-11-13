//! The title that `CWriteDB::ExtractBioseqDeflines` gives the Blast-def-line of a FASTA
//! record without `-parse_seqids`: `CFastaReader` keeps the definition line in a user
//! object as `NStr::PrintableString` writes it, and `x_GetFastaReaderDeflines` reads it
//! back (`c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp`, `c++/src/corelib/ncbistr.cpp`).

use crate::ncbi_str::parse_escapes;

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

/// The `fasta` string of `x_GetFastaReaderDeflines`: the `DefLine` field with `\t` replaced
/// by three spaces, then its escapes read (`fasta_reader_title` quotes the lines) by
/// `NStr::ParseEscapes` (`ncbi_str::parse_escapes`). The field that `CFastaReader` writes
/// (`ncbi_str::printable_string`) has no escape that fails.
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
/// of the `CFastaReader` user object (`ncbi_str::printable_string` of the definition line
/// with its `>`); `None` when the field is empty (no definition line from the user object).
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
/// No escape fails: a backslash that starts an escape is one that
/// `ncbi_str::printable_string` wrote, followed by a letter, a quote, `?`, a backslash or octal digits, or it is the
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
    use crate::ncbi_str::printable_string;

    fn title(defline: &[u8]) -> Vec<u8> {
        fasta_reader_title(&printable_string(defline))
            .unwrap()
            .unwrap()
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
