use super::*;
use crate::blastinput::ncbi_registry::{read_registry, Layer};

#[test]
fn the_protocol_and_the_browser_url_read_ncbirc() {
    assert_eq!(protocol(None), b"https:");
    assert_eq!(
        tax_browser_template(None),
        b"<@protocol@>//www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi"
    );
    let read = read_registry(
        b"[BLASTFMTUTIL]\nPROTOCOL=http:\n[TAX_BROWSER]\nBLASTFMTUTIL=<@protocol@>//x\n",
        Layer::SystemWide,
    )
    .unwrap();
    assert_eq!(protocol(Some(&read)), b"http:");
    assert_eq!(tax_browser_template(Some(&read)), b"<@protocol@>//x");
    let read = read_registry(
        b"[BLASTFMTUTIL]\nPROTOCOL=\n[TAX_BROWSER]\nBLASTFMTUTIL=\n",
        Layer::SystemWide,
    )
    .unwrap();
    // An empty PROTOCOL is set (`HasEntry` with `fCountCleared`); an empty browser URL is not
    // (`Get(...).empty()`).
    assert_eq!(protocol(Some(&read)), b"");
    assert_eq!(value(Some(&read), "TAX_BROWSER", "BLASTFMTUTIL"), None);
}

#[test]
fn a_registry_read_posts_as_ncbi() {
    // NCBI BLAST+ 2.17.0 (probe `T/t2/probe_rc`, `rep` and `bom`): the posts of one read; each
    // query's report reads `./.ncbirc` twice.
    let read = read_registry(b"[FOO]\na=1\na=2\n", Layer::SystemWide).unwrap();
    assert_eq!(
        registry_posts(&read, "blastp"),
        b"Warning: (110.8) [blastp] Found multiple [FOO] a settings; using the one from line 3\n"
            .to_vec()
    );
    let read = read_registry(b"\xef", Layer::SystemWide).unwrap();
    assert_eq!(
        registry_posts(&read, "blastp"),
        b"Error: (110.4) [blastp] Error reading the registry after line 1: \n".to_vec()
    );
}
