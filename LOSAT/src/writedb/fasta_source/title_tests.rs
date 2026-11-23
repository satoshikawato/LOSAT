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

/// NCBI's steps of `fasta_reader_title` without its search for the bytes they change, for any
/// field (also bytes that `printable_string` never writes, and escapes that fail).
fn title_by_every_step(field: &[u8]) -> Result<Option<Vec<u8>>, &'static str> {
    let fasta = fasta_reader_fasta(field)?;
    if fasta.is_empty() {
        return Ok(None);
    }
    let title = replace(&fasta[1..], b"\x01", b" ");
    Ok(Some(replace(&title, b"\t", TAB_REPLACEMENT).into_owned()))
}

#[test]
fn a_field_without_the_bytes_the_steps_change_is_its_title() {
    for text in crate::utils::byte_search::tests::search_texts() {
        let field = [b">".as_slice(), &text].concat();
        let printable = printable_string(&field);
        for field in [&text[..], &field, &printable] {
            assert_eq!(
                fasta_reader_title(field),
                title_by_every_step(field),
                "{field:?}"
            );
        }
    }
}
