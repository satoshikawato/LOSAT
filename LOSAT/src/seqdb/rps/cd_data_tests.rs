use super::*;

fn scratch(name: &str) -> PathBuf {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-cd-data-{name}-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The columns `from..to` of the profile `oid`, as a hit of one segment reads them.
fn columns(
    files: &CdProfileFiles,
    oid: usize,
    from: usize,
    to: usize,
) -> Result<CdColumns, RpsOpenError> {
    files.hit(oid, to).columns(from, to)
}

fn words(values: &[u32]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

/// A file of the profiles `data` (each a list of `Uint4`): magic, count, offsets, data.
fn file(magic: i32, profiles: &[Vec<u32>], unit: usize) -> Vec<u8> {
    let mut offsets = vec![0u32];
    for profile in profiles {
        offsets.push(offsets.last().unwrap() + (profile.len() / unit) as u32);
    }
    let mut bytes = words(&[magic as u32, profiles.len() as u32]);
    bytes.extend(words(&offsets));
    for profile in profiles {
        bytes.extend(words(profile));
    }
    bytes
}

/// Two profiles of 3 and 2 columns (and the column of zero counts that ends each): the counts of
/// row `r` of profile `p` are 1 in letter `r` and `p + 1` in letter 27.
fn write_files(dir: &Path) -> PathBuf {
    let rows = |p: u32, n: u32| -> Vec<u32> {
        let mut data = Vec::new();
        for r in 0..n {
            let mut row = vec![0u32; CD_ALPHABET];
            row[r as usize] = 1;
            row[27] = p + 1;
            data.extend(row);
        }
        // The column of zero counts that ends the profile.
        data.extend(vec![0u32; CD_ALPHABET]);
        data
    };
    let base = dir.join("cd");
    std::fs::write(
        dir.join("cd.wcounts"),
        file(RPS_MAGIC_NUM_28, &[rows(0, 3), rows(1, 2)], CD_ALPHABET),
    )
    .unwrap();
    // Profile 0: observations 2.5 for the first column, 4 for the two others (3 columns + the
    // last one, which has a run of its own); profile 1: 7 for both.
    std::fs::write(
        dir.join("cd.obsr"),
        file(
            RPS_MAGIC_NUM_28,
            &[vec![2500, 1, 4000, 3], vec![7000, 3]],
            1,
        ),
    )
    .unwrap();
    base
}

#[test]
fn the_columns_of_a_profile_are_its_counts_over_their_sums() {
    let dir = scratch("columns");
    let base = write_files(&dir);
    let files = CdProfileFiles::open(&base, "DELTABLAST").unwrap();
    let got = columns(&files, 0, 1, 3).unwrap();
    assert_eq!(got.wfreqs.len(), 2);
    // Row 1 of profile 0: one count of letter 1, one of letter 27: each 1/2.
    assert_eq!(got.wfreqs[0][1], 0.5);
    assert_eq!(got.wfreqs[0][27], 0.5);
    assert_eq!(got.wfreqs[1][2], 0.5);
    // The runs: one column of 2.5, then columns of 4.
    assert_eq!(got.iobsr, vec![4.0, 4.0]);
    let got = columns(&files, 0, 0, 2).unwrap();
    assert_eq!(got.iobsr, vec![2.5, 4.0]);
    let got = columns(&files, 1, 0, 2).unwrap();
    assert_eq!(got.wfreqs[1][1], 1.0 / 3.0);
    assert_eq!(got.wfreqs[1][27], 2.0 / 3.0);
    assert_eq!(got.iobsr, vec![7.0, 7.0]);
}

#[test]
fn a_range_past_the_profile_is_refused() {
    let dir = scratch("past");
    let base = write_files(&dir);
    let files = CdProfileFiles::open(&base, "DELTABLAST").unwrap();
    // The profile has three columns (the fourth row is the column of zero counts).
    assert!(columns(&files, 0, 0, 3).is_ok());
    assert!(matches!(
        columns(&files, 0, 0, 4),
        Err(RpsOpenError::Refused(_))
    ));
    assert!(matches!(
        columns(&files, 2, 0, 1),
        Err(RpsOpenError::Refused(_))
    ));
}

#[test]
fn the_segments_of_a_hit_read_the_columns_of_their_own_reads() {
    let dir = scratch("hit");
    let base = write_files(&dir);
    let files = CdProfileFiles::open(&base, "DELTABLAST").unwrap();
    let mut data = files.hit(0, 3);
    for (from, to) in [(1, 3), (0, 1), (2, 3), (0, 3)] {
        let hit = data.columns(from, to).unwrap();
        let own = columns(&files, 0, from, to).unwrap();
        assert_eq!(hit.wfreqs, own.wfreqs, "{from}..{to}");
        assert_eq!(hit.iobsr, own.iobsr, "{from}..{to}");
    }
    // A segment past the profile is refused as its own read is, after the others.
    assert!(matches!(data.columns(1, 4), Err(RpsOpenError::Refused(_))));
}

#[test]
fn a_hit_in_a_file_cut_short_is_refused_where_its_segment_is() {
    let dir = scratch("cut");
    let base = write_files(&dir);
    // `.wcounts` without its last two columns: profile 1 keeps its first column only.
    let path = dir.join("cd.wcounts");
    let bytes = std::fs::read(&path).unwrap();
    std::fs::write(&path, &bytes[..bytes.len() - 2 * CD_ALPHABET * 4]).unwrap();
    let files = CdProfileFiles::open(&base, "DELTABLAST").unwrap();
    let mut data = files.hit(1, 2);
    let first = data.columns(0, 1).unwrap();
    assert_eq!(first.wfreqs, columns(&files, 1, 0, 1).unwrap().wfreqs);
    let Err(RpsOpenError::Refused(text)) = data.columns(1, 2) else {
        panic!("expected LOSAT's refusal");
    };
    assert!(text.contains("is shorter than its profiles"), "{text}");
    assert!(matches!(
        columns(&files, 1, 1, 2),
        Err(RpsOpenError::Refused(_))
    ));
}

#[test]
fn a_missing_file_is_ncbis_error() {
    let dir = scratch("missing");
    let base = write_files(&dir);
    std::fs::remove_file(dir.join("cd.obsr")).unwrap();
    let Err(RpsOpenError::Ncbi(text)) = CdProfileFiles::open(&base, "DELTABLAST") else {
        panic!("expected NCBI's error");
    };
    let text = String::from_utf8_lossy(&text);
    assert!(
        text.starts_with("Cannot memory map RPS-BLAST database file: ")
            && text.ends_with("/cd.obsr"),
        "{text}"
    );
    std::fs::remove_file(dir.join("cd.wcounts")).unwrap();
    let Err(RpsOpenError::Ncbi(text)) = CdProfileFiles::open(&base, "DELTABLAST") else {
        panic!("expected NCBI's error");
    };
    let text = String::from_utf8_lossy(&text);
    assert!(text.ends_with("/cd.wcounts"), "{text}");
}

#[test]
fn a_wrong_magic_number_is_ncbis_error() {
    let dir = scratch("magic");
    let base = write_files(&dir);
    std::fs::write(dir.join("cd.wcounts"), file(7, &[vec![0; 28]], CD_ALPHABET)).unwrap();
    let Err(RpsOpenError::Ncbi(text)) = CdProfileFiles::open(&base, "DELTABLAST") else {
        panic!("expected NCBI's error");
    };
    let text = String::from_utf8_lossy(&text);
    assert!(
        text.starts_with("RPS BLAST profile file (")
            && text.ends_with(
                "/cd.wcounts) is either corrupt or constructed for an incompatible architecture"
            ),
        "{text}"
    );
}

#[test]
fn a_file_shorter_than_its_header_is_refused() {
    let dir = scratch("short");
    let base = write_files(&dir);
    std::fs::write(dir.join("cd.obsr"), [0u8; 3]).unwrap();
    assert!(matches!(
        CdProfileFiles::open(&base, "DELTABLAST"),
        Err(RpsOpenError::Refused(_))
    ));
}
