use super::*;

fn line(fmt: &str) -> Vec<u8> {
    let fields = TaxFields::parse(fmt.as_bytes()).unwrap();
    let names = TaxNames {
        taxid: 9606,
        scientific_name: b"Homo sapiens".to_vec(),
        common_name: b"human".to_vec(),
        blast_name: b"primates".to_vec(),
        super_kingdom: Ok(b"Eukaryota".to_vec()),
    };
    let mut out = Vec::new();
    fields.line(9606, &LineNames::Named(names), 3, &mut out);
    out
}

// NCBI's lines for these formats: probes T10-T12, J140-J159 of line BDC.
#[test]
fn the_format_is_read_as_ncbi_s_loop_reads_it() {
    assert_eq!(
        line("%f"),
        b"9606\tHomo sapiens\thuman\tEukaryota\tprimates\n"
    );
    assert_eq!(
        line("%f %n"),
        b"9606\tHomo sapiens\thuman\tEukaryota\tprimates\t3\n"
    );
    assert_eq!(
        line("%T %f"),
        b"9606 9606\tHomo sapiens\thuman\tEukaryota\tprimates\n"
    );
    assert_eq!(
        line("%T|%S|%L|%K|%B|%n"),
        b"9606|Homo sapiens|human|Eukaryota|primates|3\n"
    );
    assert_eq!(line("%T%n%S"), b"96063Homo sapiens\n");
    assert_eq!(line("%%T"), b"%9606\n");
    assert_eq!(line("%%x %T"), b"%%x 9606\n");
    assert_eq!(line("%x %T %n"), b"%x 9606 3\n");
    assert_eq!(line("%T %"), b"9606 %\0\n");
    assert_eq!(line("%T %%"), b"9606 %%\0\n");
    assert_eq!(line("%n"), b"3\n");
    assert!(TaxFields::parse(b"abc").is_err());
    assert!(TaxFields::parse(b"").is_err());
    assert!(TaxFields::parse(b"%x").is_err());
}

#[test]
fn missing_names_are_not_available() {
    let fields = TaxFields::parse(b"%f").unwrap();
    let mut out = Vec::new();
    fields.line(0, &LineNames::NotAvailable, 0, &mut out);
    assert_eq!(out, b"0\tN/A\tN/A\tN/A\tN/A\n");
}

// NCBI's lines for a taxid found in the index without names (probe of line RG T5b,
// `T/t5b/out/extra/tabs0.taxinfo.ncbi.out`: `9606||||`).
#[test]
fn a_taxid_found_without_names_has_empty_names() {
    let fields = TaxFields::parse(b"%T|%S|%L|%K|%B").unwrap();
    let mut out = Vec::new();
    fields.line(9606, &LineNames::Empty, 0, &mut out);
    assert_eq!(out, b"9606||||\n");
}
