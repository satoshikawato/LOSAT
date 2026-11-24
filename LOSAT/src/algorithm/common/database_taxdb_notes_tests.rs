use super::*;
use crate::report::blast_xml::HitTaxonomy;

fn lines(printer: &mut TaxdbPrinter<'_>, rank: usize, taxid: i32, hit: &HitTaxonomy) -> String {
    let mut sink = Vec::new();
    printer.at_hit((0, rank), taxid, hit, &mut sink).unwrap();
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
        placed: PlacedLoadTexts::default(),
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

/// Where the loads follow NCBI's order (a damaged taxonomy database whose last taxid shares its
/// cache slot with another taxid of a volume), a batch's prefetch texts come where the batch
/// begins (`batch_texts` before its first query's report, also where no `Hit` of the batch asks),
/// after the texts of opening, and a read's texts at the `Hit` of that query and subject, once.
#[test]
fn placed_texts_come_at_the_batch_and_the_hit_that_load() {
    let notes = TaxdbNotes {
        program: "blastp".to_string(),
        open: vec!["Error: open".to_string()],
        volumes: vec![0; 3],
        lookups: Default::default(),
        placed: PlacedLoadTexts {
            batches: vec![(0..2, b"P1\n".to_vec()), (2..3, b"P2\n".to_vec())],
            hits: [((3, 1), b"R1\n".to_vec())].into_iter().collect(),
        },
    };
    let missing = HitTaxonomy::NotFound(vec!["Error: hit".to_string()]);
    let mut printer = notes.printer();
    let at = |printer: &mut TaxdbPrinter<'_>, query: usize, rank: usize, taxid: i32| {
        let mut sink = Vec::new();
        printer
            .at_hit((query, rank), taxid, &missing, &mut sink)
            .unwrap();
        String::from_utf8(sink).unwrap()
    };
    let batch = |printer: &mut TaxdbPrinter<'_>, query: usize| {
        String::from_utf8(printer.batch_texts(query)).unwrap()
    };
    // batch 1 begins at query 0: opening, then its prefetch, before any `Hit`
    assert_eq!(batch(&mut printer, 0), "Error: [blastp] Error: open\nP1\n");
    assert_eq!(batch(&mut printer, 1), "");
    let mut sink = Vec::new();
    printer
        .at_hit((0, 2), 0, &HitTaxonomy::Unnamed, &mut sink)
        .unwrap();
    assert!(sink.is_empty());
    assert_eq!(
        at(&mut printer, 1, 0, 99990),
        "Error: [blastp] Error: hit\n"
    );
    // batch 2 begins at query 2, whose report has no `Hit` that asks
    assert_eq!(batch(&mut printer, 2), "P2\n");
    assert_eq!(batch(&mut printer, 3), "");
    assert_eq!(
        at(&mut printer, 3, 1, 99990),
        "R1\nError: [blastp] Error: hit\n"
    );
    // the read's texts are written once
    assert_eq!(
        at(&mut printer, 3, 1, 99990),
        "Error: [blastp] Error: hit\n"
    );
}
