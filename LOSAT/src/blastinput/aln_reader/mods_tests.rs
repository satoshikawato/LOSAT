use super::*;

fn ids(value: &str) -> Option<Vec<String>> {
    match seq_id_range(value.as_bytes()) {
        Ok(Some(run)) => {
            let mut ids = Vec::new();
            run.for_each(&mut |id| ids.push(String::from_utf8(id.to_vec()).unwrap()));
            Some(ids)
        }
        Ok(None) => None,
        Err(_) => panic!("{value}: refused"),
    }
}

// A secondary-accession range is spelled out with the zero padding of its first number
// (probes `f6b.m2.0xx`); a value that is not a range stays as it is.
#[test]
fn a_range_is_spelled_out() {
    assert_eq!(
        ids("AB000001-AB000003").unwrap(),
        ["AB000001", "AB000002", "AB000003"]
    );
    assert!(ids("AB000009-10").is_none());
    assert_eq!(ids("AB000001").unwrap(), ["AB000001"]);
    assert!(ids("AB0000x1").is_none());
    assert!(ids("AB000001-CD000003").is_none());
}

#[test]
fn a_backward_range_and_one_past_64_bits_are_refused() {
    assert!(matches!(
        seq_id_range(b"AB000003-AB000001"),
        Err(ModsStop::Undefined(_))
    ));
    assert!(matches!(
        seq_id_range(b"A99999999999999999999"),
        Err(ModsStop::Undefined(_))
    ));
}

// Numbers past `int` are read in 64 bits (NCBI wraps them in 32: Owner 2026-10-22), and any
// length of range is a run, not a list.
#[test]
fn numbers_past_int_and_long_ranges_are_runs() {
    assert_eq!(ids("AB9999999999").unwrap(), ["AB9999999999"]);
    assert_eq!(ids("AB2147483648").unwrap(), ["AB2147483648"]);
    assert_eq!(
        ids("AB2147483647-AB2147483649").unwrap(),
        ["AB2147483647", "AB2147483648", "AB2147483649"]
    );
    let Ok(Some(run)) = seq_id_range(b"ABCD00000001-ABCD50000000") else {
        panic!("a range of five million");
    };
    assert_eq!(run.len(), 50_000_000);
}

// The title is what the modifiers leave of the defline; an unrecognized modifier stays in it.
#[test]
fn the_title_is_the_defline_without_its_modifiers() {
    let applied = apply_defline(b"a title [organism=Homo sapiens] more", true).unwrap();
    assert_eq!(applied.title.as_deref(), Some(&b"a title more"[..]));
    assert!(applied.effects.content.has_any());
    let applied = apply_defline(b"[organism=Homo sapiens]", true).unwrap();
    assert!(applied.title.is_none());
}

/// No modifier value can panic: every value that a table of this file gives to a descriptor is
/// a value of the ENUMERATED type it is written as (NCBI's `fli_cdna` is `fli-cdna` in
/// `seq.asn`, which a spelling in this file once broke), and the descriptors of each are made.
#[test]
fn the_tables_name_values_of_their_enumerations() {
    let mut deflines: Vec<String> = Vec::new();
    let text = |key: &[u8]| String::from_utf8_lossy(key).into_owned();
    for name in SUBSOURCE.iter().chain(ORGMOD) {
        deflines.push(format!("[{}=x]", text(name)));
        deflines.push(format!("[{}=true]", text(name)));
    }
    for (key, _, _) in BIOMOL {
        deflines.push(format!("[mol-type={}]", text(key)));
    }
    for (key, _) in TECH {
        deflines.push(format!("[tech={}]", text(key)));
    }
    for (key, _) in COMPLETENESS {
        deflines.push(format!("[completeness={}]", text(key)));
    }
    for origin in [
        "natural",
        "natural mutant",
        "mutant",
        "artificial",
        "synthetic",
        "?",
    ] {
        deflines.push(format!("[origin={origin}]"));
    }
    for location in [
        "genomic",
        "chloroplast",
        "chromoplast",
        "kinetoplast",
        "mitochondrial",
        "mitochondrion",
        "plastid",
        "macronuclear",
        "extrachromosomal",
        "plasmid",
        "transposon",
        "insertion sequence",
        "cyanelle",
        "proviral",
        "provirus",
        "virion",
        "nucleomorph",
        "apicoplast",
        "leucoplast",
        "proplastid",
        "endogenous virus",
        "hydrogenosome",
        "chromosome",
        "chromatophore",
        "plasmid in mitochondrion",
        "plasmid in plastid",
        "unknown",
    ] {
        deflines.push(format!("[location={location}]"));
    }
    for defline in &deflines {
        for protein in [true, false] {
            if let Ok(applied) = apply_defline(defline.as_bytes(), protein) {
                let _ = applied.effects.content.descriptors();
            }
        }
    }
    // Each table entry is accepted (the defline is not refused and writes a descriptor), so the
    // loop above made a value for it.
    for (key, name) in TECH {
        let applied = apply_defline(format!("[tech={}]", text(key)).as_bytes(), true).unwrap();
        assert!(applied.effects.content.has_any(), "tech {name}");
    }
    for (key, name) in COMPLETENESS {
        let applied =
            apply_defline(format!("[completeness={}]", text(key)).as_bytes(), true).unwrap();
        assert!(applied.effects.content.has_any(), "completeness {name}");
    }
    for (key, name, _) in BIOMOL {
        let applied = apply_defline(format!("[mol-type={}]", text(key)).as_bytes(), false).unwrap();
        assert!(applied.effects.content.has_any(), "mol-type {name}");
    }
    for name in SUBSOURCE {
        assert!(subsource_kind(name).is_ok(), "{}", text(name));
    }
    for name in ORGMOD {
        assert!(orgmod_kind(name).is_ok(), "{}", text(name));
    }
    for (_, name, _) in BIOMOL {
        assert!(checked(named::MOLINFO, "biomol", name).is_ok(), "{name}");
    }
    for (_, name) in TECH {
        assert!(checked(named::MOLINFO, "tech", name).is_ok(), "{name}");
    }
    for (_, name) in COMPLETENESS {
        assert!(
            checked(named::MOLINFO, "completeness", name).is_ok(),
            "{name}"
        );
    }
    for (name, member) in [
        ("natmut", "origin"),
        ("mut", "origin"),
        ("mitochondrion", "genome"),
        ("proviral", "genome"),
        ("extrachrom", "genome"),
        ("insertion-seq", "genome"),
    ] {
        assert!(checked(named::BIOSOURCE, member, name).is_ok(), "{name}");
    }
}
