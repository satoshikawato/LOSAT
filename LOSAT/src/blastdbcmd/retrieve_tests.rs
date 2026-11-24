use super::*;

#[test]
fn version5_keys_are_preprocessed_as_ncbi_s() {
    let key = preprocess_accession_for_v5;
    assert_eq!(key(b"p69905"), b"P69905");
    assert_eq!(key(b"np_000002.1"), b"NP_000002.1");
    assert_eq!(key(b"sp|P99999.2|TEST_HUMAN"), b"P99999.2");
    assert_eq!(key(b"pdb|1ABC|A"), b"1ABC_A");
    assert_eq!(key(b"28zx_m"), b"28zx_m");
    assert_eq!(key(b"lcl|myprot9"), b"myprot9");
    assert_eq!(key(b"gnl|mydb|prot10"), b"mydb:prot10");
    assert_eq!(key(b"gnl|BL_ORD_ID|3"), b"BL_ORD_ID:3");
    assert_eq!(key(b"pir||s00007"), b"pir||s00007");
    assert_eq!(key(b"gi|123456"), b"123456");
    assert_eq!(key(b"mydb:prot10"), b"MYDB:PROT10");
}

/// The answers of `x_GetOids` in the fixture database `name` (version 4) of
/// `tests/blastdbcmd/db`.
fn v4_oids(name: &str) -> impl Fn(&str) -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/blastdbcmd/db")
        .join(name);
    let db = SeqDb::open(path.to_str().unwrap(), SeqType::Protein).unwrap();
    assert!(!db.is_version5());
    move |id: &str| {
        let ids = DbIds {
            lookup: Some(db.id_lookup().unwrap()),
            version5: false,
        };
        match get_oids(&ids, id.as_bytes()) {
            Lookup::Oids(oids) => format!("oids {oids:?}"),
            Lookup::NotFound(acc) => format!("not found {}", String::from_utf8_lossy(&acc)),
            Lookup::Stop(Stop::Ncbi(message)) => format!("error {message}"),
            Lookup::Stop(Stop::NcbiBytes(message)) => {
                format!("error {}", String::from_utf8_lossy(&message))
            }
            Lookup::Stop(stop) => format!("{stop:?}"),
        }
    }
}

/// NCBI's version check of an accession (`x_CheckVersions`) reads the header of each
/// candidate (`GetSeqIDs`), and its read of a record with a UTF8String title posts NCBI's
/// warning before the lookup's own message (review 2 of RG-T5d, D2): `IdLookup::take_utf8`
/// names the read; a lookup of the stored version reads no header.
#[cfg(unix)]
#[test]
fn the_version_check_of_an_accession_notes_its_utf8_string_reads() {
    let source = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastdbcmd/db");
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-vcheck-utf8-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    for ext in ["pin", "psq", "psd", "psi"] {
        std::fs::copy(
            source.join(format!("d4.{ext}")),
            dir.join(format!("d4.{ext}")),
        )
        .unwrap();
    }
    // The title tag of OID 1 (`XP_800001.1 dup one a`) made a UTF8String's.
    let index = std::fs::read(dir.join("d4.pin")).unwrap();
    let word = |at: usize| u32::from_be_bytes(index[at..at + 4].try_into().unwrap()) as usize;
    let mut at = 8;
    for _ in 0..2 {
        at += 4 + word(at);
    }
    let first_header = at + 4 + 8 + 4;
    let mut header = std::fs::read(source.join("d4.phr")).unwrap();
    let title = word(first_header + 4) + 6;
    assert_eq!(header[title], 0x1a);
    header[title] = 0x0c;
    std::fs::write(dir.join("d4.phr"), header).unwrap();
    let db = SeqDb::open(dir.join("d4").to_str().unwrap(), SeqType::Protein).unwrap();
    let lookup = db.id_lookup().unwrap();
    // A version that no sequence holds: the candidates' headers (`XP_800001.1` and `.2`, OIDs 1
    // and 256) are read, the first with the warning.
    assert!(accession_to_oids_v4(&lookup, b"XP_800001.3")
        .unwrap()
        .is_empty());
    let note = lookup.take_utf8().expect("the read of OID 1");
    assert_eq!(
        db.take_utf8_post_at(note.at).as_deref(),
        Some("Blast-def-line-set.E.title")
    );
    assert!(lookup.take_utf8().is_none());
    // The stored version: no header is read.
    assert_eq!(
        accession_to_oids_v4(&lookup, b"XP_800001.1").unwrap(),
        vec![1]
    );
    assert!(lookup.take_utf8().is_none());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn ordinal_ids_of_version4_are_read_without_files() {
    let oids = v4_oids("p4n");
    assert_eq!(oids("gnl|BL_ORD_ID|3"), "oids [3]");
    assert_eq!(oids("gnl|BL_ORD_ID|3|"), "oids [3]");
    assert_eq!(oids("gnl|BL_ORD_ID|3\r"), "oids [3]");
    assert_eq!(oids("gnl|BL_ORD_ID|999"), "not found gnl|BL_ORD_ID|999");
    assert_eq!(oids(" gnl|BL_ORD_ID|1"), "not found  gnl|BL_ORD_ID|1");
    assert_eq!(oids("GNL|bl_ord_id|3"), "not found GNL|bl_ord_id|3");
    assert_eq!(oids("3"), "not found 3");
    assert_eq!(oids("myprot9"), "not found myprot9");
    assert_eq!(oids(""), "not found ");
    for tag in ["0", "-1", "abc", "03"] {
        assert_eq!(
            oids(&format!("gnl|BL_ORD_ID|{tag}")),
            "error CObject_id::GetId(): Invalid choice selection: NCBI-General::Object-id.str"
        );
    }
}

#[test]
fn ids_of_a_parsed_version4_database_are_looked_up_in_its_isams() {
    let oids = v4_oids("p4");
    assert_eq!(oids("gnl|BL_ORD_ID|3"), "oids [3]");
    assert_eq!(oids("XP_000020.1"), oids("xp_000020.1"));
    assert!(oids("XP_000020.1").starts_with("oids ["));
    assert_eq!(oids("1819"), oids("gi|1819"));
    assert!(oids("1819").starts_with("oids ["));
    assert_eq!(oids("NOSUCH9"), "not found NOSUCH9");
}

#[test]
fn batch_lines_split_in_two() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat_bdc_batch_{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("b.txt");
    std::fs::write(
        &file,
        "gnl|BL_ORD_ID|3\n  gnl|BL_ORD_ID|4\t3-\t\n\n   \nx 1-5 21\r\n",
    )
    .unwrap();
    let lines = read_batch(
        "entry_batch",
        file.as_os_str(),
        &mut Cin::new(&mut std::io::empty()),
    )
    .unwrap();
    // The `\r` of `\r\n` stays on Linux and goes on Windows and macOS (`NcbiGetlineEOL`).
    let crlf_rest = match LineEnds::HOST {
        LineEnds::Unix => "1-5 21\r",
        LineEnds::Windows | LineEnds::Darwin => "1-5 21",
    };
    assert_eq!(
        lines,
        [
            ("gnl|BL_ORD_ID|3", ""),
            ("gnl|BL_ORD_ID|4", "3-"),
            ("x", crlf_rest),
        ]
        .map(|(id, rest)| (id.as_bytes().to_vec(), rest.as_bytes().to_vec()))
    );
    std::fs::remove_dir_all(&dir).unwrap();
    assert!(
        matches!(read_batch("entry_batch", std::ffi::OsStr::new("in/nosuch.txt"), &mut Cin::new(&mut std::io::empty())),
            Err(Stop::Ncbi(message))
        if message == "Argument \"entry_batch\". File is not accessible:  `in/nosuch.txt'")
    );
}
