//! The title that `CWriteDB::ExtractBioseqDeflines` gives the Blast-def-line of a FASTA
//! record without `-parse_seqids`: `CFastaReader` keeps the definition line in a user
//! object as `NStr::PrintableString` writes it, and `x_GetFastaReaderDeflines` reads it
//! back (`c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp`, `c++/src/corelib/ncbistr.cpp`).

use crate::ncbi_str::parse_escapes;

/// `NStr::Replace` and `NStr::ReplaceInPlace`: every occurrence of `search`, from the
/// left, without reading the replacement again.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3371-3385
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
const TAB_REPLACEMENT: &[u8] = b"   ";

/// The title of `x_GetFastaReaderDeflines` without `parse_ids`, from the `DefLine` field
/// of the `CFastaReader` user object (`ncbi_str::printable_string` of the definition line
/// with its `>`); `None` when the field is empty (no definition line from the user object).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1752-1760,1766-1767,1782-1798
/// The field is searched for the two characters `\t` before its escapes are read, so a
/// tab of the definition line (written `\t`) becomes three spaces, and so does a
/// backslash followed by `t` (written `\\t`): the search finds the second backslash with
/// the `t`, and the first backslash then escapes the first of the three spaces. A ctrl-A
/// becomes one space. Every other byte comes back as it was read.
///
/// No escape fails: a backslash that starts an escape is one that `ncbi_str::printable_string`
/// wrote, followed by a letter, a quote, `?`, a backslash or octal digits, or it is the backslash
/// left of a replaced `\\t`, followed by a space.
pub(super) fn fasta_reader_title(defline_field: &[u8]) -> Result<Option<Vec<u8>>, &'static str> {
    let fasta = fasta_reader_fasta(defline_field)?;
    if fasta.is_empty() {
        return Ok(None);
    }
    let title = replace(&fasta[1..], b"\x01", b" ");
    Ok(Some(replace(&title, b"\t", TAB_REPLACEMENT)))
}

#[cfg(test)]
#[path = "title_tests.rs"]
mod tests;
