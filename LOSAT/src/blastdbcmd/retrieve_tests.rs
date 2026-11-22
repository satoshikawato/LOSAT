use super::*;

#[test]
fn version5_keys_are_preprocessed_as_ncbi_s() {
    let key = preprocess_accession_for_v5;
    assert_eq!(key("p69905"), "P69905");
    assert_eq!(key("np_000002.1"), "NP_000002.1");
    assert_eq!(key("sp|P99999.2|TEST_HUMAN"), "P99999.2");
    assert_eq!(key("pdb|1ABC|A"), "1ABC_A");
    assert_eq!(key("28zx_m"), "28zx_m");
    assert_eq!(key("lcl|myprot9"), "myprot9");
    assert_eq!(key("gnl|mydb|prot10"), "mydb:prot10");
    assert_eq!(key("gnl|BL_ORD_ID|3"), "BL_ORD_ID:3");
    assert_eq!(key("pir||s00007"), "pir||s00007");
    assert_eq!(key("gi|123456"), "123456");
    assert_eq!(key("mydb:prot10"), "MYDB:PROT10");
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
        match get_oids(&ids, id) {
            Lookup::Oids(oids) => format!("oids {oids:?}"),
            Lookup::NotFound(acc) => format!("not found {acc}"),
            Lookup::Stop(Stop::Ncbi(message)) => format!("error {message}"),
            Lookup::Stop(Stop::NcbiBytes(message)) => {
                format!("error {}", String::from_utf8_lossy(&message))
            }
            Lookup::Stop(stop) => format!("{stop:?}"),
        }
    }
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
        file.to_str().unwrap(),
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
            ("gnl|BL_ORD_ID|3".to_string(), String::new()),
            ("gnl|BL_ORD_ID|4".to_string(), "3-".to_string()),
            ("x".to_string(), crlf_rest.to_string()),
        ]
    );
    std::fs::remove_dir_all(&dir).unwrap();
    assert!(
        matches!(read_batch("entry_batch", "in/nosuch.txt", &mut Cin::new(&mut std::io::empty())),
            Err(Stop::Ncbi(message))
        if message == "Argument \"entry_batch\". File is not accessible:  `in/nosuch.txt'")
    );
}
