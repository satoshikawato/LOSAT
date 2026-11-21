use super::*;

fn tokens(text: &str) -> Vec<String> {
    let mut warnings = Vec::new();
    read_memory_si_list(text.as_bytes(), &mut warnings)
        .into_iter()
        .map(|t| String::from_utf8(t).unwrap())
        .collect()
}

#[test]
fn text_tokens_as_ncbi_reads_them() {
    assert_eq!(
        tokens(">AAA00001.1\n# comment line\naaa00001.2 ref|NP_000001.1| # trailing\r\nlcl|x\tgnl|db|t\n"),
        ["AAA00001.1", "aaa00001.2", "ref|NP_000001.1|", "lcl|x", "gnl|db|t"]
    );
    assert_eq!(tokens("a#1 b"), ["a#1", "b"]);
    assert_eq!(tokens(">>  >id\r\n"), ["id"]);
    assert_eq!(tokens("P69905"), ["P69905"]);
    assert!(tokens("# only\n\n").is_empty());
    assert!(tokens("").is_empty());
}

#[test]
fn a_token_of_vertical_tabs_is_reported_and_dropped() {
    let mut warnings = Vec::new();
    let ids = read_memory_si_list(b"\x0b\x0c ok\n", &mut warnings);
    assert_eq!(ids, [b"ok".to_vec()]);
    assert_eq!(
        warnings,
        ["WARNING:  \x0b\x0c is not a valid seqid string."]
    );
}

fn binary(ids: &[&[u8]], db_vol_length: u64) -> Vec<u8> {
    let mut body = Vec::new();
    body.extend_from_slice(&(ids.len() as u64).to_le_bytes());
    body.extend_from_slice(&5u32.to_le_bytes());
    body.extend_from_slice(b"title");
    body.push(3);
    body.extend_from_slice(b"Oct");
    body.extend_from_slice(&db_vol_length.to_le_bytes());
    if db_vol_length != 0 {
        body.push(2);
        body.extend_from_slice(b"Dt");
        body.extend_from_slice(&4u32.to_le_bytes());
        body.extend_from_slice(b"vols");
    }
    for id in ids {
        if id.len() >= 0xFF {
            body.push(0xFF);
            body.extend_from_slice(&(id.len() as u32).to_le_bytes());
        } else {
            body.push(id.len() as u8);
        }
        body.extend_from_slice(id);
    }
    let size = 1 + 8 + body.len();
    let mut file = vec![0u8];
    file.extend_from_slice(&(size as u64).to_le_bytes());
    file.extend_from_slice(&body);
    file
}

#[test]
fn a_binary_list_gives_its_ids_and_header() {
    let long = vec![b'A'; 300];
    let file = binary(&[b"AAA00001.2", b"Q00002", &long], 12345);
    let mut warnings = Vec::new();
    let list = read_si_list(&file, &mut warnings).unwrap();
    assert_eq!(
        list.ids,
        [b"AAA00001.2".to_vec(), b"Q00002".to_vec(), long.clone()]
    );
    assert_eq!(
        list.info,
        SeqidListInfo {
            is_v4: false,
            file_size: file.len() as u64,
            num_ids: 3,
            title: b"title".to_vec(),
            create_date: b"Oct".to_vec(),
            db_vol_length: 12345,
            db_create_date: b"Dt".to_vec(),
            db_vol_names: b"vols".to_vec(),
        }
    );
    assert!(warnings.is_empty());
    let plain = binary(&[b"X"], 0);
    let list = read_si_list(&plain, &mut warnings).unwrap();
    assert_eq!(list.ids, [b"X".to_vec()]);
    assert_eq!(list.info.db_vol_length, 0);
    assert!(list.info.db_vol_names.is_empty());
}

#[test]
fn damaged_binary_lists_give_ncbis_messages() {
    let mut warnings = Vec::new();
    assert_eq!(
        read_si_list(b"", &mut warnings).unwrap_err(),
        "Failed to map seqidlist file "
    );
    let mut file = binary(&[b"X"], 0);
    file.push(0);
    assert_eq!(
        read_si_list(&file, &mut warnings).unwrap_err(),
        "Invalid seqidlist file"
    );
    let mut short = binary(&[b"X", b"Y"], 0);
    let cut = short.len() - 2;
    short.truncate(cut);
    short[1..9].copy_from_slice(&(cut as u64).to_le_bytes());
    assert_eq!(
        read_si_list(&short, &mut warnings).unwrap_err(),
        "Invalid total num of ids in seqidlist file"
    );
    // A binary list of no ids is read as text: its header bytes are the tokens.
    let empty = binary(&[], 0);
    let list = read_si_list(&empty, &mut warnings).unwrap();
    assert!(list.info.is_v4);
    assert!(!list.ids.is_empty());
}

fn processed_warned(tokens: &[&str]) -> (Vec<String>, Vec<SeqIdWarning>) {
    let ids: Vec<Vec<u8>> = tokens.iter().map(|t| t.as_bytes().to_vec()).collect();
    let mut warnings = Vec::new();
    let ids = process_for_v5(&ids, &mut |warning| warnings.push(warning))
        .unwrap()
        .into_iter()
        .map(|t| String::from_utf8(t).unwrap())
        .collect();
    (ids, warnings)
}

fn processed(tokens: &[&str]) -> Vec<String> {
    let (ids, warnings) = processed_warned(tokens);
    assert_eq!(warnings, []);
    ids
}

#[test]
fn version5_ids_are_rewritten_as_ncbi_rewrites_them() {
    assert_eq!(
        processed(&[
            "aaa00001.2",
            "ref|NP_000001.1|",
            "np_000002",
            "lcl|myprot9",
            "gnl|mydb|tag1",
            "mydb:tag1",
            "12345",
            "gi|999",
            "sp|Q00001.1|N1_HUMAN",
            "MixedCase1",
            "pdb|1ABC|A",
            "1ABC_A",
            "1ABC",
            "pir||S00007",
            "AAA00001.2",
        ]),
        [
            "12345",
            "1ABC",
            "1ABC_A",
            "AAA00001.2",
            "MixedCase1",
            "NP_000001.1",
            "NP_000002",
            "Q00001.1",
            "mydb:tag1",
            "myprot9",
            "pir||S00007",
        ]
    );
    assert!(processed(&["gi|999", "gi|998"]).is_empty());
    // An unknown FASTA tag is a local id under `fParse_AnyLocal` (probe S23 of
    // dropin-survey-20261029/DI8: no error, the id is looked up as it is).
    assert_eq!(processed(&["xyz|ABC|DEF"]), ["xyz|ABC|DEF"]);
    let bad = process_for_v5(&[b"gi|abc".to_vec()], &mut |_| {}).unwrap_err();
    assert_eq!(bad.message(), b"Invalid seq id: gi|abc");
    // A FASTA-style id of several parts: NCBI's warning for each non-empty part after the
    // first id (rows db.gisl.bp.full.6, db.gisl.bn.full.6), the gi dropped.
    let warning = SeqIdWarning {
        subcode: 10,
        message: b"Ignoring extra parts (synonyms?) in FASTA-style ID gi|1111|ref|NP_1.1|".to_vec(),
    };
    assert_eq!(
        processed_warned(&["gi|1111|ref|NP_1.1|", "NP_2.1"]),
        (vec!["NP_2.1".to_string()], vec![warning.clone(), warning])
    );
}
