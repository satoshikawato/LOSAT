use super::*;

fn found(name: &str, seqtype: SeqType) -> FoundDb {
    FoundDb {
        quoted: format!("\"{name}\"").into_bytes(),
        seqtype,
    }
}

// Probes L01-L05, J160-J178 of line BDC.
#[test]
fn the_format_is_read_as_ncbi_s_loop_reads_it() {
    let format = ListFormat::parse(b"%f %p").unwrap();
    assert_eq!(format.offsets, [0, 3]);
    assert_eq!(format.types, b"fp");
    let format = ListFormat::parse(b"%%f %f %%").unwrap();
    assert_eq!(format.spec, b"%f %f %");
    assert_eq!(format.offsets, [3]);
    assert_eq!(format.replace(&[b"db/x".to_vec()]), b"%f db/x %");
    let format = ListFormat::parse(b"%f %").unwrap();
    assert_eq!(format.types, b"f\0");
    assert_eq!(format.replace(&[b"a".to_vec(), b"b".to_vec()]), b"a b");
    // After `%%` the loop steps over the byte that follows the second `%`.
    assert!(ListFormat::parse(b"%%%f").is_err());
    assert!(ListFormat::parse(b"abc").is_err());
    assert!(ListFormat::parse(b"").is_err());
    assert!(ListFormat::parse(b"%%").is_err());
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

#[test]
fn the_extensions_follow_the_version() {
    assert_eq!(file_extensions(false).len(), 20);
    let v5 = file_extensions(true);
    assert_eq!(v5.len(), 25);
    assert!(!v5.contains(&"si"));
    assert!(v5.contains(&"tf-lock"));
}
