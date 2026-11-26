use super::*;
use crate::report::blast_xml::HitTaxonomy;

fn lines(printer: &mut TaxdbPrinter<'_>, taxid: i32, hit: &HitTaxonomy) -> String {
    let mut sink = Vec::new();
    printer.at_hit(taxid, hit, &mut sink).unwrap();
    String::from_utf8(sink).unwrap()
}

/// NCBI writes the texts of opening that a `Hit`'s lookup made at the first `Hit` with a taxid
/// (a negative one too: the volume's lookup of a load asks for positive taxids only, the `Hit`
/// for every taxid that is not 0), and the `Hit`'s own lookup every time, each after the lines
/// before it; the loads' texts are the prefetch's or the reads' (`read_reported_taxonomy`).
#[test]
fn the_texts_of_opening_come_at_the_first_hit_that_asks() {
    let notes = TaxdbNotes {
        program: "blastp".to_string(),
        open: vec!["Error: open".to_string()],
    };
    let mut printer = notes.printer();
    let missing = HitTaxonomy::NotFound(vec!["Error: hit".to_string()]);
    // a subject without a taxid writes nothing and does not open the database
    assert_eq!(lines(&mut printer, 0, &HitTaxonomy::Unnamed), "");
    assert_eq!(
        lines(&mut printer, -7, &missing),
        "Error: [blastp] Error: open\nError: [blastp] Error: hit\n"
    );
    // the database is open: the `Hit`'s own lookup alone
    assert_eq!(
        lines(&mut printer, 9606, &missing),
        "Error: [blastp] Error: hit\n"
    );
    assert_eq!(
        lines(&mut printer, 9606, &HitTaxonomy::Named(Vec::new())),
        ""
    );
    // no texts (a load opened the database, or no taxonomy database)
    let none = TaxdbNotes::default();
    assert_eq!(lines(&mut none.printer(), 9606, &HitTaxonomy::Unnamed), "");
}
