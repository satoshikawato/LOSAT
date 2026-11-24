use super::*;

fn found(name: &str, seqtype: SeqType) -> FoundDb {
    FoundDb {
        quoted: format!("\"{name}\"").into_bytes(),
        seqtype,
    }
}

#[test]
fn the_list_sorts_on_the_quoted_names() {
    let mut dbs = [
        found("d/a!", SeqType::Nucleotide),
        found("d/zz", SeqType::Protein),
        found("d/x.00", SeqType::Protein),
        found("d/a b", SeqType::Protein),
        found("d/x", SeqType::Protein),
        found("d/zz", SeqType::Nucleotide),
        found("d/.h", SeqType::Protein),
    ];
    let rank = |seqtype: SeqType| matches!(seqtype, SeqType::Nucleotide) as u8;
    dbs.sort_by(|a, b| {
        a.quoted
            .cmp(&b.quoted)
            .then(rank(a.seqtype).cmp(&rank(b.seqtype)))
    });
    let names: Vec<String> = dbs
        .iter()
        .map(|db| String::from_utf8(db.file_name()).unwrap())
        .collect();
    assert_eq!(
        names,
        ["d/.h", "d/a b", "d/a!", "d/x", "d/x.00", "d/zz", "d/zz"]
    );
    assert_eq!(dbs[5].seqtype, SeqType::Protein);
}

#[test]
fn a_dblist_name_removes_every_name_that_holds_it() {
    assert!(contains(b"\"db/p5\"", b"\"db/p5\""));
    assert!(!contains(b"\"db/p5n\"", b"\"db/p5\""));
    assert!(!contains(b"\"db/j/lst/p5\"", b"\"p5\""));
    assert!(contains(b"\"db/n55\"", b"n5"));
    assert!(contains(b"\"db/s18mv.00\"", b"s18mv.00"));
}
