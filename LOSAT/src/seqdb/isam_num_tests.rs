use super::*;

/// Writes the files of `rows` (sorted) in the given key width; the index samples every
/// `page`-th row.
fn write_isam(dir: &Path, rows: &[(i64, i32)], long_id: bool, page: usize) -> std::path::PathBuf {
    let base = dir.join("v");
    let mut data = Vec::new();
    let mut samples = Vec::new();
    for (at, &(key, oid)) in rows.iter().enumerate() {
        let mut row = Vec::new();
        if long_id {
            row.extend_from_slice(&key.to_be_bytes());
        } else {
            row.extend_from_slice(&(key as u32).to_be_bytes());
        }
        row.extend_from_slice(&oid.to_be_bytes());
        if at % page == 0 {
            samples.extend_from_slice(&row);
        }
        data.extend_from_slice(&row);
    }
    let term = if long_id { 12 } else { 8 };
    let mut index = Vec::new();
    for word in [
        ISAM_VERSION,
        if long_id { NUMERIC_LONG_ID } else { NUMERIC },
        data.len() as i32,
        rows.len() as i32,
        (samples.len() / term) as i32,
        page as i32,
        0,
        0,
        0,
    ] {
        index.extend_from_slice(&word.to_be_bytes());
    }
    index.extend_from_slice(&samples);
    index.extend_from_slice(&vec![0xff; term - 4]);
    index.extend_from_slice(&0i32.to_be_bytes());
    std::fs::write(with_extension(&base, SeqType::Protein, "ni"), index).unwrap();
    std::fs::write(with_extension(&base, SeqType::Protein, "nd"), data).unwrap();
    base
}

fn temp_dir(stem: &str) -> std::path::PathBuf {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat_isam_num_{stem}_{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn every_row_is_found_in_both_key_widths() {
    for long_id in [false, true] {
        let dir = temp_dir(if long_id { "long" } else { "short" });
        let step = if long_id { 1i64 << 33 } else { 7 };
        let rows: Vec<(i64, i32)> = (1..=23).map(|i| (i * step, i as i32 + 100)).collect();
        let base = write_isam(&dir, &rows, long_id, 4);
        let isam = NumericIsam::open_volume(&base, SeqType::Protein)
            .unwrap()
            .unwrap();
        for &(key, oid) in &rows {
            assert_eq!(isam.id_to_oid(key).unwrap(), Some(oid), "key {key}");
            assert_eq!(isam.id_to_oid(key + 1).unwrap(), None, "key {}", key + 1);
            assert_eq!(isam.id_to_oid(key - 1).unwrap(), None, "key {}", key - 1);
        }
        assert_eq!(isam.id_to_oid(-5).unwrap(), None);
        assert_eq!(isam.id_to_oid(i64::MAX).unwrap(), None);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}

#[test]
fn missing_files_are_silent_and_a_bad_header_is_refused() {
    let dir = temp_dir("bad");
    let base = write_isam(&dir, &[(5, 1), (9, 2)], false, 256);
    let nd = with_extension(&base, SeqType::Protein, "nd");
    let ni = with_extension(&base, SeqType::Protein, "ni");
    let data = std::fs::read(&nd).unwrap();
    std::fs::write(&nd, &data[..8]).unwrap();
    assert!(NumericIsam::open_volume(&base, SeqType::Protein)
        .unwrap()
        .is_none());
    std::fs::remove_file(&nd).unwrap();
    assert!(NumericIsam::open_volume(&base, SeqType::Protein)
        .unwrap()
        .is_none());
    std::fs::write(&nd, &data).unwrap();
    let mut index = std::fs::read(&ni).unwrap();
    index[3] = 2;
    std::fs::write(&ni, &index).unwrap();
    assert!(matches!(
        NumericIsam::open_volume(&base, SeqType::Protein),
        Err(error @ SeqDbError::NcbiCrashes { .. }) if error.ncbi_crashes()
            && error.to_string().contains(": unsupported database: a GI index file of version 2 and type ")
            && error.to_string().ends_with(" (NCBI crashes on it)")
    ));
    std::fs::write(&ni, &index[..20]).unwrap();
    assert!(matches!(
        NumericIsam::open_volume(&base, SeqType::Protein),
        Err(SeqDbError::NcbiCrashes { .. })
    ));
    std::fs::remove_dir_all(&dir).unwrap();
}

fn list(keys: &[i64]) -> Vec<NumOid> {
    keys.iter().map(|&key| NumOid { key, oid: None }).collect()
}

fn opened(base: &Path) -> NumericIsam {
    match NumericIsam::open_for_list(base, SeqType::Protein, ListIsam::Gi, 1).unwrap() {
        Some(Opened::Ready(isam)) => isam,
        other => panic!("not ready: {other:?}"),
    }
}

#[test]
fn a_list_is_translated_against_every_page_of_the_volume() {
    for long_id in [false, true] {
        let dir = temp_dir(if long_id { "batch_long" } else { "batch_short" });
        let rows: Vec<(i64, i32)> = (1..=23).map(|i| (i * 7, i as i32)).collect();
        let base = write_isam(&dir, &rows, long_id, 4);
        let isam = opened(&base);
        let mut ids = list(&[1, 7, 8, 28, 28, 70, 161, 162]);
        isam.ids_to_oids(100, &mut ids).unwrap();
        let oids: Vec<Option<u32>> = ids.iter().map(|id| id.oid).collect();
        assert_eq!(
            oids,
            vec![
                None,
                Some(101),
                None,
                Some(104),
                Some(104),
                Some(110),
                Some(123),
                None
            ]
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}

#[test]
fn an_id_translated_by_an_earlier_volume_keeps_its_oid() {
    let dir = temp_dir("batch_kept");
    let base = write_isam(&dir, &[(5, 0), (9, 1), (12, 2)], false, 2);
    let isam = opened(&base);
    let mut ids = list(&[5, 9, 12]);
    ids[1].oid = Some(3);
    isam.ids_to_oids(10, &mut ids).unwrap();
    let oids: Vec<Option<u32>> = ids.iter().map(|id| id.oid).collect();
    assert_eq!(oids, vec![Some(10), Some(3), Some(12)]);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_list_isam_needs_oids_and_an_index_of_the_data_length() {
    let dir = temp_dir("batch_open");
    let base = write_isam(&dir, &[(5, 0), (9, 1)], false, 256);
    assert!(
        NumericIsam::open_for_list(&base, SeqType::Protein, ListIsam::Gi, 0)
            .unwrap()
            .is_none()
    );
    assert!(
        NumericIsam::open_for_list(&base, SeqType::Protein, ListIsam::Ti, 2)
            .unwrap()
            .is_none()
    );
    let nd = with_extension(&base, SeqType::Protein, "nd");
    let data = std::fs::read(&nd).unwrap();
    std::fs::write(&nd, &data[..8]).unwrap();
    assert!(matches!(
        NumericIsam::open_for_list(&base, SeqType::Protein, ListIsam::Gi, 2).unwrap(),
        Some(Opened::Uninitialized)
    ));
    std::fs::remove_dir_all(&dir).unwrap();
}
