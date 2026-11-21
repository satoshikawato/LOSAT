use super::*;
use crate::report::blast_xml::HitTaxonomy;

fn lines(printer: &mut TaxdbPrinter<'_>, rank: usize, taxid: i32, hit: &HitTaxonomy) -> String {
    let mut sink = Vec::new();
    printer.at_hit(rank, taxid, hit, &mut sink).unwrap();
    String::from_utf8(sink).unwrap()
}

/// NCBI writes the texts of opening at the first `Hit` with a taxid, the volume's lookup
/// at the first `Hit` of each (volume, taxid), and the `Hit`'s own lookup every time, each
/// after the lines before it (and so after the reader's messages of the batch).
#[test]
fn the_notes_come_at_the_first_hit_that_needs_them() {
    let notes = TaxdbNotes {
        program: "blastp".to_string(),
        open: vec!["Error: open".to_string()],
        volumes: vec![0, 0, 1, 0],
        lookups: [
            ((0, 9606), vec!["Error: lookup 9606".to_string()]),
            ((1, 9606), vec!["Error: lookup 9606".to_string()]),
        ]
        .into_iter()
        .collect(),
    };
    let mut printer = notes.printer();
    let missing = HitTaxonomy::NotFound(vec!["Error: hit".to_string()]);
    // a subject without a taxid writes nothing and does not open the database
    assert_eq!(lines(&mut printer, 3, 0, &HitTaxonomy::Unnamed), "");
    assert_eq!(
        lines(&mut printer, 0, 9606, &missing),
        "Error: [blastp] Error: open\nError: [blastp] Error: lookup 9606\nError: [blastp] Error: hit\n"
    );
    // the same (volume, taxid): the volume does not look it up again, the `Hit` does
    assert_eq!(
        lines(&mut printer, 1, 9606, &missing),
        "Error: [blastp] Error: hit\n"
    );
    assert_eq!(
        lines(&mut printer, 1, 9606, &HitTaxonomy::Named(Vec::new())),
        ""
    );
    // another volume: its own lookup
    assert_eq!(
        lines(&mut printer, 2, 9606, &HitTaxonomy::Named(Vec::new())),
        "Error: [blastp] Error: lookup 9606\n"
    );
    // no notes (a batch after the first, or a database without a taxonomy database)
    let none = TaxdbNotes::default();
    assert_eq!(
        lines(&mut none.printer(), 0, 9606, &HitTaxonomy::Unnamed),
        ""
    );
}
