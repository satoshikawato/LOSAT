//! The BlastOutput documents of the fixtures (`tests/fixtures/rpsbproc/arch/x_*.xml.gz`, NCBI
//! rpsblast and rpstblastn 2.17.0 `-outfmt 5`; NCBI's rpsbproc 0.5.1 output in `oracle/`):
//! the queries that `ParseBlastOutput` makes of them.

use super::*;
use crate::rpsbproc::data::DataFiles;

fn fixture(path: &str) -> String {
    format!(
        "{}/tests/fixtures/rpsbproc/{path}",
        env!("CARGO_MANIFEST_DIR")
    )
}

/// The data files of `dir` (`data` or `data_with/<name>`, which holds all six).
fn data(dir: &str) -> CddData {
    let file = |name: &str| fixture(&format!("{dir}/{name}"));
    CddData::load(&DataFiles {
        cddids: file("cddid.tbl"),
        spfeats: file("cddannot.dat"),
        genfeats: file("cddannot_generic.dat"),
        bscores: file("bitscore_specific.txt"),
        cdtrack: file("cdtrack.txt"),
        famlinks: file("family_superfamily_links"),
    })
    .unwrap()
}

/// The documents of `arch/<name>.xml.gz` with rpsbproc's default cutoff.
fn documents(name: &str, data: &CddData) -> Vec<Document> {
    let file = std::fs::File::open(fixture(&format!("arch/{name}.xml.gz"))).unwrap();
    let mut stream = XmlStream::new(flate2::read::GzDecoder::new(file), b"");
    let mut documents = Vec::new();
    while let Some(document) = stream.next_document(0.01, data).unwrap() {
        documents.push(document);
    }
    documents
}

fn query(document: &Document, index: usize) -> &DomQuery {
    document.queries[index].as_ref().unwrap()
}

/// Whether the hit of `pssm_id` is a specific hit of frame `frame` (`IsSpecific`, among the
/// frame's concise or full set).
fn specific(query: &DomQuery, frame: usize, pssm_id: u32) -> bool {
    let index = query
        .aligns
        .iter()
        .position(|align| align.pssm_id == pssm_id)
        .unwrap();
    let align = &query.aligns[index];
    query.frames[frame].sorted.contains(&index) && align.rep && align.spec_qualified
}

// `x_special`: NCBI's reader decodes the entity references of the titles, also the `&amp;` and
// `&lt;x&gt;` that NCBI's writer left as text (oracle `x.special.both.rep`: `sp2 & <x> ]]> -->
// <!-- x`).
#[test]
fn titles_with_entity_references() {
    let data = data("data_with/xnew");
    let documents = documents("x_special", &data);
    assert_eq!(documents.len(), 1);
    let document = &documents[0];
    assert_eq!(
        (&document.version[..], &document.db[..]),
        (&b"RPSBLAST 2.17.0+"[..], &b"Cdd"[..])
    );
    assert_eq!(document.matrix.as_deref(), Some(&b"BLOSUM62"[..]));
    assert_eq!(document.expect, 0.01);
    assert_eq!(query(document, 0).title, b"sp1 a&b <c> \"d\" 'e'".to_vec());
    assert_eq!(
        query(document, 1).title,
        b"sp2 & <x> ]]> --> <!-- x".to_vec()
    );
}

// `x_ntspec1` and `x_ntspec2`: hits whose segments come from the alignment strings
// (`ParseAlignStrings`) are specific where the same hits of the archive path, from its
// `Std-seg`s, are not (PS16d probe 1, `ana2full.txt`): PRK00915.1-1 (PSSM-ID 503205) on frame -2
// of mg00017 (oracle `x.ntspec1.both.rep`) and COG1047 (440668) on frame +1 of mg00006 with
// `-lcase_masking` (oracle `x.ntspec2.both.full`).
#[test]
fn hits_specific_by_their_alignment_strings() {
    let data = data("data_with/xnew");
    let spec1 = documents("x_ntspec1", &data);
    let query1 = query(&spec1[0], 0);
    assert!(query1.is_na && query1.translated.is_some());
    assert!(specific(query1, 4, 503205));
    let spec2 = documents("x_ntspec2", &data);
    assert!(specific(query(&spec2[0], 0), 0, 440668));
}

// `x_cattrunc` (`x_one`, then the first half of `x_three`): the cut document does not read and
// ends the documents; `x_ntnohit`: a nucleotide query without hits is a protein query (no
// HSP commits it; oracle `x.ntnohit.both.rep`: `Peptide`, no `DOMAINS`).
#[test]
fn documents_until_one_does_not_read() {
    let data = data("data");
    assert_eq!(documents("x_cattrunc", &data).len(), 1);
    assert_eq!(documents("x_cat2", &data).len(), 2);
    let nohit = documents("x_ntnohit", &data_for_nt());
    let query = query(&nohit[0], 0);
    assert!(!query.is_na && query.aligns.is_empty() && query.translated.is_none());
}

/// The data of the nucleotide fixtures (`data_with/nt`).
fn data_for_nt() -> CddData {
    data("data_with/nt")
}
