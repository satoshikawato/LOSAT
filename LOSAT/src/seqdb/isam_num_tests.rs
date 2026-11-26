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

fn temp_dir(stem: &str) -> crate::utils::test_dir::TempDir {
    crate::utils::test_dir::TempDir::unique(&format!("isam-num-{stem}"))
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
}

/// NCBI's constructor reads the first and last rows of the data file whatever `x_InitSearch`
/// returned (`x_FindIndexBounds`): inside the last memory page of the mapping it reads the zeros
/// after the file, past it memory that is not the file's, which LOSAT refuses in both states.
#[test]
fn the_bounds_are_read_in_the_mapping_of_the_data_file() {
    let dir = temp_dir("bounds");
    let base = write_isam(&dir, &[(5, 0), (9, 1)], false, 256);
    let ni = with_extension(&base, SeqType::Protein, "ni");
    let index = std::fs::read(&ni).unwrap();
    let with_terms = |terms: i32, data_len: i32| {
        let mut index = index.clone();
        index[8..12].copy_from_slice(&data_len.to_be_bytes());
        index[12..16].copy_from_slice(&terms.to_be_bytes());
        std::fs::write(&ni, index).unwrap();
    };
    let open = || NumericIsam::open_for_list(&base, SeqType::Protein, ListIsam::Gi, 2);
    // 512 rows of 8 bytes end at the end of the first page: zeros past the 16 bytes of the file.
    with_terms(512, 16);
    let Some(Opened::Ready(isam)) = open().unwrap() else {
        panic!("not ready");
    };
    assert_eq!((isam.first_key, isam.last_key), (5, 0));
    with_terms(512, 4096);
    assert!(matches!(open().unwrap(), Some(Opened::Uninitialized)));
    // One more row ends past the page: refused whether or not the lengths agree.
    for data_len in [16, 4104] {
        with_terms(513, data_len);
        let error = open().unwrap_err();
        assert!(error.ncbi_crashes(), "{error}");
        assert!(
            error
                .to_string()
                .contains("NCBI reads row 512 of the GI data file (its last of 513 terms), past the memory page after its 16 bytes"),
            "{error}"
        );
    }
    // A header with more samples than the rows fill: the last page would start past the last
    // row, and the last key is still the row before `m_NumTerms` (review 2 of line RG D1c-1, R3).
    let base3 = write_isam(&dir, &[(5, 0), (9, 1), (12, 2)], false, 2);
    let mut more_samples = std::fs::read(&ni).unwrap();
    more_samples[16..20].copy_from_slice(&3i32.to_be_bytes());
    std::fs::write(&ni, more_samples).unwrap();
    let Some(Opened::Ready(isam)) =
        NumericIsam::open_for_list(&base3, SeqType::Protein, ListIsam::Gi, 3).unwrap()
    else {
        panic!("not ready");
    };
    assert_eq!((isam.first_key, isam.last_key), (5, 12));
    // An empty data file has no mapping (NCBI's null pointer): its first row is refused too.
    std::fs::write(with_extension(&base, SeqType::Protein, "nd"), []).unwrap();
    with_terms(1, 0);
    assert!(open().unwrap_err().ncbi_crashes());
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
}

#[test]
fn a_negative_scan_marks_excluded_oids_visible_and_unlisted_oids_included() {
    let dir = temp_dir("negative_scan");
    // OID 0 has GIs 10 and 20; OID 1 has GI 30; OID 2 has GI 40 (outside vol_end=2? no, vol_end=10).
    let rows = [(10, 0), (20, 0), (30, 1), (40, 2), (50, 3)];
    let base = write_isam(&dir, &rows, false, 2);
    let isam = opened(&base);
    let mut out = super::super::isam::NegativeOids::new();
    // Exclude GI 10 (OID 0 still has unlisted GI 20 -> stays included!) and GI 30 (OID 1 only has GI 30 -> excluded!).
    isam.negative_scan(10, 14, &[10, 30, 99], &mut out).unwrap();
    assert!(out.get_oid_status(10));
    assert!(!out.get_oid_status(11));
    assert!(out.get_oid_status(12));
    assert!(out.get_oid_status(13));
}

#[test]
fn list_isam_pig_opens_ppi_and_ppd_when_both_exist() {
    let dir = temp_dir("pig_open");
    let base = write_isam(&dir, &[(1001, 0), (1002, 1)], false, 256);
    assert!(
        NumericIsam::open_for_list(&base, SeqType::Protein, ListIsam::Pig, 2)
            .unwrap()
            .is_none()
    );
    std::fs::copy(
        with_extension(&base, SeqType::Protein, "ni"),
        with_extension(&base, SeqType::Protein, "pi"),
    )
    .unwrap();
    std::fs::copy(
        with_extension(&base, SeqType::Protein, "nd"),
        with_extension(&base, SeqType::Protein, "pd"),
    )
    .unwrap();
    let Some(Opened::Ready(isam)) =
        NumericIsam::open_for_list(&base, SeqType::Protein, ListIsam::Pig, 2).unwrap()
    else {
        panic!("PIG ISAM not ready");
    };
    let mut ids = list(&[1001, 1002, 9999]);
    isam.ids_to_oids(5, &mut ids).unwrap();
    assert_eq!(
        ids.iter().map(|id| id.oid).collect::<Vec<_>>(),
        vec![Some(5), Some(6), None]
    );
}
