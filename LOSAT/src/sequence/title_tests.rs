use super::*;

fn no_source(ids: &[ParsedSeqId], title: &[u8]) -> Vec<u8> {
    get_title_without_source(ids, title)
}

fn named(ids: &[ParsedSeqId], title: &[u8], name: Option<&str>) -> Vec<u8> {
    let name = name.map(|name| Arc::<[u8]>::from(name.as_bytes()));
    let Ok(text) = get_title(ids, title, true, || Ok::<_, std::convert::Infallible>(name));
    text
}

#[test]
fn titles_lose_trailing_dots_and_spaces_and_third_party_ids_get_tpa() {
    let tpg = ParsedSeqId::Tpg(Default::default());
    assert_eq!(no_source(&[], b"rec4. ."), b"rec4");
    assert_eq!(
        no_source(std::slice::from_ref(&tpg), b"third party"),
        b"TPA: third party"
    );
    assert_eq!(no_source(std::slice::from_ref(&tpg), b"tpa: x"), b"tpa: x");
    assert_eq!(no_source(&[tpg], b""), b"");
}

#[test]
fn titles_keep_html_entities_and_other_punctuation() {
    // `CDeflineGenerator` would decode the entities and trim `,;~`: GetTitle does neither.
    assert_eq!(no_source(&[], b"a&amp;b"), b"a&amp;b");
    assert_eq!(no_source(&[], b"abc,"), b"abc,");
    assert_eq!(no_source(&[], b"semi;"), b"semi;");
    assert_eq!(no_source(&[], b" lead"), b" lead");
    assert_eq!(no_source(&[], b"..."), b"");
}

#[test]
fn an_empty_title_with_a_source_is_its_organism_or_no_definition_line() {
    assert_eq!(named(&[], b"", Some("homo sapiens")), b"Homo sapiens");
    assert_eq!(named(&[], b". .", Some(" \tE. coli \n")), b"E. coli");
    assert_eq!(named(&[], b"", None), b"No definition line found");
    assert_eq!(named(&[], b"", Some(" ")), b"No definition line found");
    assert_eq!(named(&[], b"titled", None), b"titled");
    let tpg = ParsedSeqId::Tpg(Default::default());
    assert_eq!(named(&[tpg], b"", None), b"TPA: No definition line found");
}

#[test]
fn the_name_is_asked_only_for_a_title_made_from_the_source() {
    let asked = std::cell::Cell::new(false);
    let ask = || {
        asked.set(true);
        Err::<Option<Arc<[u8]>>, ()>(())
    };
    assert_eq!(get_title(&[], b"titled", true, ask), Ok(b"titled".to_vec()));
    assert!(!asked.get());
    assert_eq!(get_title(&[], b"", false, ask), Ok(Vec::new()));
    assert!(!asked.get());
    assert_eq!(get_title(&[], b"", true, ask), Err(()));
    assert!(asked.get());
}
