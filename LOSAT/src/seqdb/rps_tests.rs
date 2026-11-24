use super::*;

/// A temporary directory of this test (tests of one binary run in parallel).
fn scratch(name: &str) -> PathBuf {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-rps-{name}-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn words(values: &[i32]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

/// A volume of two profiles (lengths 3 and 2), with `columns` columns, whose table holds
/// the word ending at row 2 in a cell and three words in one overflow run.
pub(crate) fn write_volume(dir: &Path, columns: usize, aux: &str) -> PathBuf {
    let base = dir.join("v");
    let magic = if columns == 26 {
        RPS_MAGIC_NUM
    } else {
        RPS_MAGIC_NUM_28
    };
    let offsets = [0, 4, 7];
    let mut rps = vec![magic, 2];
    rps.extend_from_slice(&offsets);
    for row in 0..7 {
        for column in 0..columns {
            rps.push((row * 100 + column) as i32);
        }
    }
    std::fs::write(base.with_extension("rps"), words(&rps)).unwrap();
    let mut loo = vec![magic, 1, 0, 0, 0, 0, 0, 0, 40, 0];
    let mut cells = vec![[0i32; 4]; RPS_BACKBONE_SIZE + 1];
    cells[5] = [1, 2, 0, 0];
    cells[9] = [4, 2, 4, 0];
    for cell in &cells {
        loo.extend_from_slice(cell);
    }
    loo.extend_from_slice(&[99, 3, 5, 6]);
    std::fs::write(base.with_extension("loo"), words(&loo)).unwrap();
    std::fs::write(base.with_extension("aux"), aux).unwrap();
    base
}

const AUX: &str =
    "BLOSUM62\n11\n1\n0.000000e+00\n0.000000e+00\n0\n0\n100.000000\n3\n4.1e-02\n2\n3.9e-02\n";

#[test]
fn the_aux_text_is_read_as_the_streams_do() {
    let aux = parse_aux(AUX).unwrap();
    assert_eq!(aux.matrix, "BLOSUM62");
    assert_eq!((aux.gap_open, aux.gap_extend, aux.scale), (11, 1, 100.0));
    assert_eq!(aux.karlin_k, vec![0.041, 0.039]);
    // Any white space separates the tokens.
    assert_eq!(parse_aux(&AUX.replace('\n', " \t")).unwrap(), aux);
    // A pair without its K, and a token of the wrong kind.
    assert!(parse_aux(&format!("{AUX}7\n")).is_none());
    assert!(parse_aux(&AUX.replace("\n11\n", "\neleven\n")).is_none());
}

#[test]
fn a_volume_of_28_columns_is_viewed_in_place() {
    let dir = scratch("v28");
    let base = write_volume(&dir, 28, AUX);
    let files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
    assert_eq!(files.alphabet, 28);
    assert_eq!(files.offsets(), &[0, 4, 7]);
    assert_eq!(files.rows().len(), 7);
    assert_eq!(files.rows()[5][27], 527);
    assert_eq!(files.cells().unwrap()[9], [4, 2, 4, 0]);
    assert_eq!(files.overflow().unwrap(), &[99, 3, 5, 6]);
    assert_eq!(files.pv()[0], (1 << 5) | (1 << 9));
    let _ = std::fs::remove_dir_all(&dir);
}

// Review C of line RG D1a-4 (SF2 c): NCBI's row after the rows (`rps_pssm[num_pssm_rows]`,
// `row_after_rows`) is the file's own bytes after the rows, then the zeros of its last memory
// page, and a letter whose bytes lie past that page is `None` (NCBI's read faults: SIGBUS).
#[test]
fn the_row_after_the_rows_reads_the_files_bytes_then_its_last_pages_zeros() {
    let dir = scratch("past");
    let base = write_volume(&dir, 28, AUX);
    let path = base.with_extension("rps");
    let rows = std::fs::read(&path).unwrap();
    // The file ends with the rows (byte 804): the page's zeros.
    let files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
    assert_eq!(files.row_after_rows(), [Some(0); 28]);
    drop(files);
    // A row's 112 bytes after the rows are read as written; 8 bytes: two letters, then zeros.
    let after = words(&(1000..1028).collect::<Vec<i32>>());
    for (bytes, read) in [(112, 28), (8, 2)] {
        std::fs::write(&path, [rows.as_slice(), &after[..bytes]].concat()).unwrap();
        let files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
        let expected: [Option<i32>; 28] = std::array::from_fn(|letter| {
            Some(if letter < read {
                1000 + letter as i32
            } else {
                0
            })
        });
        assert_eq!(
            files.row_after_rows(),
            expected,
            "{bytes} bytes after the rows"
        );
        drop(files);
    }
    // Three profiles of 36 rows, whose rows and file end 40 bytes before the end of the page
    // (24 + 36 * 112 = 4056 bytes): letters 0 to 9 read its zeros, 10 to 27 lie past it.
    let mut page = vec![RPS_MAGIC_NUM_28, 3, 0, 12, 24, 36];
    page.extend((0..36 * 28).map(|value: i32| value + 1));
    std::fs::write(&path, words(&page)).unwrap();
    assert_eq!(std::fs::metadata(&path).unwrap().len(), 4056);
    std::fs::write(
        base.with_extension("aux"),
        AUX.replace(
            "3\n4.1e-02\n2\n3.9e-02\n",
            "11\n4.1e-02\n11\n3.9e-02\n11\n4.0e-02\n",
        ),
    )
    .unwrap();
    let files = RpsVolumeFiles::open(&base, &[11, 11, 11]).unwrap();
    let row = files.row_after_rows();
    assert_eq!(row[..10], [Some(0); 10]);
    assert_eq!(row[10..], [None; 18]);
    drop(files);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_volume_of_26_columns_reads_the_next_rows_first_columns() {
    let dir = scratch("v26");
    let base = write_volume(&dir, 26, AUX);
    let files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
    assert_eq!(files.alphabet, 26);
    let rows = files.rows();
    assert_eq!(rows[2][25], 225);
    // Letters 26 and 27 read the next row's columns 0 and 1 (stride 26).
    assert_eq!((rows[2][26], rows[2][27]), (300, 301));
    assert_eq!((rows[6][26], rows[6][27]), (0, 0));
    let _ = std::fs::remove_dir_all(&dir);
}

/// A `.freq` file of `columns` columns for the two profiles of `write_volume` (rows 0 to
/// 6): the ratio of row `r`, column `c` is `r * 1000 + c` over 10^9.
fn write_freq(base: &Path, columns: usize) {
    let magic = if columns == 26 {
        RPS_MAGIC_NUM
    } else {
        RPS_MAGIC_NUM_28
    };
    let mut freq = vec![magic, 2, 0, 4, 7];
    for row in 0..7 {
        for column in 0..columns {
            freq.push((row * 1000 + column) as i32);
        }
    }
    std::fs::write(base.with_extension("freq"), words(&freq)).unwrap();
}

#[test]
fn the_freq_rows_have_the_rps_alphabets_stride() {
    let dir = scratch("freq");
    let mut ratios = Vec::new();
    let mut bytes = Vec::new();
    // 28 columns on both sides: rows as written, ÷ 10^9.
    let base = write_volume(&dir, 28, AUX);
    write_freq(&base, 28);
    let mut files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
    files.open_freq(&base, "RPSBLAST", None).unwrap();
    assert!(files.freq_ratios(4, 2, &mut bytes, &mut ratios).unwrap());
    assert_eq!(ratios.len(), 2);
    assert_eq!(ratios[0][0], 4000.0 / 1e9);
    assert_eq!(ratios[1][27], 5027.0 / 1e9);
    // The next volume is written over these files. (Windows does not write a mapped file.)
    drop(files);
    // A 26-column `.rps` with a 28-column `.freq` reads rows two integers early: row r
    // starts at integer 26 r of the rows (NCBI's drift), and columns 26 and 27 are 0.0.
    let base = write_volume(&dir, 26, AUX);
    write_freq(&base, 28);
    let mut files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
    files.open_freq(&base, "RPSBLAST", None).unwrap();
    assert!(files.freq_ratios(2, 1, &mut bytes, &mut ratios).unwrap());
    let flat: Vec<i32> = (0..7)
        .flat_map(|row| (0..28).map(move |column| row * 1000 + column))
        .collect();
    assert_eq!(ratios[0][0], f64::from(flat[52]) / 1e9);
    assert_eq!(ratios[0][25], f64::from(flat[77]) / 1e9);
    assert_eq!((ratios[0][26], ratios[0][27]), (0.0, 0.0));
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn the_score_bound_of_a_profile_sums_its_rows_positive_maxima() {
    let dir = scratch("bounds");
    let base = write_volume(&dir, 28, AUX);
    let files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
    // Row r holds r * 100 + c: the maximum of a row is r * 100 + 27. Profile 0 has rows
    // 0..=3 (the last is its sentinel row), profile 1 rows 4..=6.
    let bounds = files.profile_score_bounds(&base).unwrap();
    assert_eq!(bounds, vec![27 + 127 + 227 + 327, 427 + 527 + 627]);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_missing_or_bad_freq_file_is_ncbis_error_and_a_short_one_is_refused() {
    let dir = scratch("freqerr");
    let base = write_volume(&dir, 28, AUX);
    let shown = canonical_volume_path(&base);
    let text = format!(
        "{shown} contains no frequency ratios needed for composition-based statistics.\nPlease disable composition-based statistics when searching against {shown}."
    );
    let mut files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
    assert_eq!(
        files.open_freq(&base, "RPSBLAST", None).err(),
        Some(RpsOpenError::Ncbi(text.clone()))
    );
    std::fs::write(base.with_extension("freq"), words(&[7, 2, 0, 4, 7])).unwrap();
    assert_eq!(
        files.open_freq(&base, "RPSBLAST", None).err(),
        Some(RpsOpenError::Ncbi(text))
    );
    // A database of one volume: its `-db` value as given (`RpsVolumeSetup::with_database`).
    assert_eq!(
        files.open_freq(&base, "RPSTBLASTN", Some("./cdd")).err(),
        Some(RpsOpenError::Ncbi(
            "./cdd contains no frequency ratios needed for composition-based statistics.\nPlease disable composition-based statistics when searching against ./cdd.".to_string()
        ))
    );
    write_freq(&base, 28);
    let freq = std::fs::read(base.with_extension("freq")).unwrap();
    std::fs::write(base.with_extension("freq"), &freq[..freq.len() - 4]).unwrap();
    assert!(matches!(
        files.open_freq(&base, "RPSBLAST", None),
        Err(RpsOpenError::Refused(_))
    ));
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn ncbis_errors_come_in_ncbis_file_order() {
    let dir = scratch("errors");
    let base = write_volume(&dir, 28, AUX);
    let shown = canonical_volume_path(&base);
    let open = || RpsVolumeFiles::open(&base, &[3, 2]).err();
    std::fs::write(base.with_extension("rps"), words(&[7, 2])).unwrap();
    assert_eq!(
        open(),
        Some(RpsOpenError::Ncbi(format!(
            "RPS BLAST profile file ({shown}.rps) is either corrupt or constructed for an incompatible architecture"
        )))
    );
    std::fs::remove_file(base.with_extension("loo")).unwrap();
    assert_eq!(
        open(),
        Some(RpsOpenError::Ncbi(format!(
            "Cannot memory map RPS-BLAST database file: {shown}.loo"
        )))
    );
    std::fs::remove_file(base.with_extension("aux")).unwrap();
    assert_eq!(
        open(),
        Some(RpsOpenError::Ncbi(format!(
            "Cannot open RPS-BLAST auxiliary file: {shown}.aux"
        )))
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn inconsistent_files_are_refused() {
    let dir = scratch("refused");
    let base = write_volume(&dir, 28, AUX);
    // `.pin` lengths other than the profiles' columns open: NCBI reads them only in the traceback,
    // where `rps_search::profile` asks each read (series PR, item PS10e).
    let files = RpsVolumeFiles::open(&base, &[3, 3]).unwrap();
    drop(files);
    // A profile without its sentinel row: offsets that do not ascend.
    let rps = std::fs::read(base.with_extension("rps")).unwrap();
    let mut flat = rps.clone();
    flat[12..16].copy_from_slice(&0i32.to_le_bytes());
    std::fs::write(base.with_extension("rps"), &flat).unwrap();
    let Err(RpsOpenError::Refused(text)) = RpsVolumeFiles::open(&base, &[3, 2]) else {
        panic!("offsets that do not ascend are refused");
    };
    assert!(
        text.contains(
            ".rps profile 0 has no rows, not even its sentinel row (its offsets 0 and 0 do not ascend)"
        ),
        "{text}"
    );
    // Offsets that descend (a hand-edited `.rps`): the same text, no row count below 0.
    let mut descending = rps.clone();
    descending[12..16].copy_from_slice(&(-3i32).to_le_bytes());
    std::fs::write(base.with_extension("rps"), &descending).unwrap();
    let Err(RpsOpenError::Refused(text)) = RpsVolumeFiles::open(&base, &[3, 2]) else {
        panic!("offsets that descend are refused");
    };
    assert!(
        text.contains(
            ".rps profile 0 has no rows, not even its sentinel row (its offsets 0 and -3 do not ascend)"
        ),
        "{text}"
    );
    std::fs::write(base.with_extension("rps"), &rps).unwrap();
    // Too few K values.
    std::fs::write(base.with_extension("aux"), AUX.replace("2\n3.9e-02\n", "")).unwrap();
    assert!(matches!(
        RpsVolumeFiles::open(&base, &[3, 2]),
        Err(RpsOpenError::Refused(_))
    ));
    // A .rps cut inside its rows (NCBI crashes).
    let base = write_volume(&dir, 28, AUX);
    let rps = std::fs::read(base.with_extension("rps")).unwrap();
    std::fs::write(base.with_extension("rps"), &rps[..rps.len() - 4]).unwrap();
    assert!(matches!(
        RpsVolumeFiles::open(&base, &[3, 2]),
        Err(RpsOpenError::Refused(_))
    ));
    // Bytes after the last row are not read (NCBI searches the file).
    let mut longer = rps.clone();
    longer.extend_from_slice(&[0x7f; 112]);
    std::fs::write(base.with_extension("rps"), &longer).unwrap();
    let files = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
    assert_eq!(files.rows().len(), 7);
    // (Windows does not write a mapped file.)
    drop(files);
    // An empty .rps (NCBI crashes).
    let base = write_volume(&dir, 28, AUX);
    std::fs::write(base.with_extension("rps"), b"").unwrap();
    assert!(matches!(
        RpsVolumeFiles::open(&base, &[3, 2]),
        Err(RpsOpenError::Refused(_))
    ));
    let _ = std::fs::remove_dir_all(&dir);
}

// D6 of PS12s: the files of a volume checked once in the run open again without LOSAT's check of
// the table's hits (`open_checked`): the presence bits are the same, and a hit outside the
// rows, which `open_for` refuses, is not read again.
#[test]
fn files_checked_in_the_run_open_without_the_check_of_the_hits() {
    let dir = scratch("checked");
    let base = write_volume(&dir, 28, AUX);
    let checked = RpsVolumeFiles::open(&base, &[3, 2]).unwrap();
    let reopened = RpsVolumeFiles::open_checked(&base, &[3, 2], "RPSBLAST", false).unwrap();
    assert_eq!(reopened.pv(), checked.pv());
    assert_eq!(reopened.rows(), checked.rows());
    drop((checked, reopened));
    // The hit of word 5 at row 99 of 7 (cell 5's first entry, byte 40 + 5 * 16 + 4).
    let mut loo = std::fs::read(base.with_extension("loo")).unwrap();
    loo[124..128].copy_from_slice(&99i32.to_le_bytes());
    std::fs::write(base.with_extension("loo"), &loo).unwrap();
    assert!(matches!(
        RpsVolumeFiles::open(&base, &[3, 2]),
        Err(RpsOpenError::Refused(text)) if text.contains(".loo word 5 names a position")
    ));
    let reopened = RpsVolumeFiles::open_checked(&base, &[3, 2], "RPSBLAST", false).unwrap();
    assert_ne!(reopened.pv()[0] & (1 << 5), 0);
    drop(reopened);
    // NCBI's own errors of the files stay: a missing .loo.
    std::fs::remove_file(base.with_extension("loo")).unwrap();
    assert!(matches!(
        RpsVolumeFiles::open_checked(&base, &[3, 2], "RPSBLAST", false),
        Err(RpsVolumeFault::Files(RpsOpenError::Ncbi(_)))
    ));
    let _ = std::fs::remove_dir_all(&dir);
}

/// A volume of two profiles of 140,000 and 100 letters (140,102 rows: 69 buckets, so two groups
/// of `RPS_GROUP_BUCKETS`) whose table holds the given cells' hits (word, hits in the cell's
/// order), the long lists in the overflow array.
fn write_wide_volume(dir: &Path, cell_hits: &[(usize, &[i32])]) -> PathBuf {
    let base = dir.join("wide");
    let (a, b) = (140_000usize, 100usize);
    let mut rps = vec![RPS_MAGIC_NUM_28, 2, 0, (a + 1) as i32, (a + b + 2) as i32];
    rps.extend(std::iter::repeat_n(3i32, (a + b + 2) * 28));
    std::fs::write(base.with_extension("rps"), words(&rps)).unwrap();
    let mut cells = vec![[0i32; 4]; RPS_BACKBONE_SIZE + 1];
    let mut overflow = Vec::new();
    for &(word, hits) in cell_hits {
        let mut cell = [hits.len() as i32, 0, 0, 0];
        if hits.len() as i32 <= RPS_HITS_PER_CELL {
            cell[1..=hits.len()].copy_from_slice(hits);
        } else {
            cell[1] = hits[0];
            cell[2] = (overflow.len() * 4) as i32;
            overflow.extend_from_slice(&hits[1..]);
        }
        cells[word] = cell;
    }
    let mut loo = vec![RPS_MAGIC_NUM_28, 1, 0, 0, 0, 0, 0, 0, 40, 0];
    loo.extend(cells.iter().flatten());
    loo.extend(overflow);
    std::fs::write(base.with_extension("loo"), words(&loo)).unwrap();
    std::fs::write(base.with_extension("aux"), AUX).unwrap();
    base
}

// PS12m: the walk of `check_cells` (`hits_in_rows_and_groups`) and the regroup's test
// (`list_is_grouped`) are one predicate: they agree on random cells, ordered and not, in
// a cell's head and run, with entries in and out of the rows.
#[test]
fn the_check_walk_and_the_regroup_test_agree_on_grouped_cells() {
    let mut state = 12345u64;
    let mut next = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as u32
    };
    let (mut grouped, mut ungrouped) = (0, 0);
    for case in 0..4_000 {
        let len = 1 + next() as usize % 40;
        let rows = 140_000u32;
        let mut list: Vec<i32> = (0..len).map(|_| 2 + (next() % (rows - 2)) as i32).collect();
        if case % 3 != 0 {
            // Ascending by group (what makeprofiledb writes), the order inside a group free.
            list.sort_by_key(|&entry| entry_group(entry, RPS_GROUP_BUCKETS));
        }
        if case % 11 == 0 {
            list[len / 2] = 1; // a hit outside the rows
        }
        let split = 1.min(len);
        let (head, run) = list.split_at(split);
        let mut walk_ungrouped = false;
        let in_rows = hits_in_rows_and_groups(head, run, rows - 2, &mut walk_ungrouped);
        assert_eq!(in_rows, list.iter().all(|&e| (2..rows as i32).contains(&e)));
        assert_eq!(
            !walk_ungrouped,
            list_is_grouped(head, run, RPS_GROUP_BUCKETS),
            "{list:?}"
        );
        if walk_ungrouped {
            ungrouped += 1;
        } else {
            grouped += 1;
        }
    }
    assert!(grouped > 500 && ungrouped > 500, "{grouped} {ungrouped}");
}

// PS12m: the groups of a cell's hits (entry - 2, per 2,048 rows, per 64 buckets) are
// non-decreasing in every table that makeprofiledb writes; a table that lists a cell otherwise
// is regrouped by a stable sort of the cell (the hits of one group keep their order), the copy
// read in place of the file's cells and overflow array; the other cells, and the cells whose
// buckets are unordered inside a group, are not touched, and grouped files are not copied.
#[test]
fn cells_not_in_groups_of_buckets_are_regrouped_stably() {
    let dir = scratch("regroup");
    const G1: i32 = 2 + 65 * 2048; // the first entry of bucket 65: group 1
                                   // Group 1 is entries from 2 + 64 * 2048 = 131,074.
    let drop_in_cell = [G1 + 5, 10];
    let drop_after_first = [G1 + 5, 10, 20, G1 + 7]; // the drop is between the cell and the run
    let drop_in_run = [10, G1, 20, G1 + 9, 30];
    let grouped = [10, 20, G1];
    let unordered_in_group = [3000, 10, 2500, 40_000, G1, G1 + 1]; // buckets drop, groups do not
    let cells: [(usize, &[i32]); 5] = [
        (5, &drop_in_cell),
        (9, &drop_in_run),
        (11, &grouped),
        (13, &unordered_in_group),
        (15, &drop_after_first),
    ];
    let base = write_wide_volume(&dir, &cells);
    let lengths = [140_000, 100];
    let files = RpsVolumeFiles::open(&base, &lengths).unwrap();
    assert!(files.is_regrouped());
    let (cells, overflow) = (files.cells().unwrap(), files.overflow().unwrap());
    assert_eq!(cells[5], [2, 10, G1 + 5, 0]);
    // The long list: 10 20 30 (group 0), then G1 and G1 + 9 in their order.
    assert_eq!(&cells[9][..2], &[5, 10]);
    let (head, run) = cell_hits(&cells[9], overflow).unwrap();
    assert_eq!((head, run), (&[10][..], &[20, 30, G1, G1 + 9][..]));
    assert_eq!(cells[11], [3, 10, 20, G1]);
    // Unordered buckets inside group 0 stay as they are.
    let (head, run) = cell_hits(&cells[13], overflow).unwrap();
    assert_eq!(
        (head, run),
        (&[3000][..], &[10, 2500, 40_000, G1, G1 + 1][..])
    );
    // The drop between the cell's first hit and the run's first.
    let (head, run) = cell_hits(&cells[15], overflow).unwrap();
    assert_eq!((head, run), (&[10][..], &[20, G1 + 5, G1 + 7][..]));
    // The reopening of a file checked in the run: regrouped when told so, as it was.
    let reopened = RpsVolumeFiles::open_checked(&base, &lengths, "RPSBLAST", true).unwrap();
    assert_eq!(reopened.cells(), files.cells());
    assert_eq!(reopened.overflow(), files.overflow());
    let raw = RpsVolumeFiles::open_checked(&base, &lengths, "RPSBLAST", false).unwrap();
    assert!(!raw.is_regrouped());
    assert_eq!(raw.cells().unwrap()[5], [2, G1 + 5, 10, 0]);
    drop((files, reopened, raw));
    // A table whose cells are grouped (a bucket drop inside a group is not a drop) is read in place.
    let cells: [(usize, &[i32]); 2] = [
        (9, &[10, 20, 3000, 2500, 40_000]),
        (13, &unordered_in_group),
    ];
    let base = write_wide_volume(&dir, &cells);
    let files = RpsVolumeFiles::open(&base, &lengths).unwrap();
    assert!(!files.is_regrouped());
    assert_eq!(files.cells().unwrap()[9][1], 10);
    let _ = std::fs::remove_dir_all(&dir);
}

// The sort of one cell: stable, so a group's hits keep the order they had, and idempotent on
// a grouped list; the regrouping works at any group size.
#[test]
fn a_cells_hits_are_sorted_by_group_and_keep_their_order_inside_one() {
    // Buckets (entry - 2) / 2048: 2, 0, 1, 0, 2, 1, 0.
    let original = [5000, 12, 2100, 30, 5001, 2050, 20];
    for (group_buckets, sorted) in [
        (1, [12, 30, 20, 2100, 2050, 5000, 5001]),
        (2, [12, 2100, 30, 2050, 20, 5000, 5001]),
    ] {
        let mut hits = original;
        assert!(!list_is_grouped(&hits[..1], &hits[1..], group_buckets));
        regroup_cell(&mut hits, group_buckets);
        assert_eq!(hits, sorted, "{group_buckets} buckets a group");
        assert!(list_is_grouped(&hits[..1], &hits[1..], group_buckets));
        let again = hits;
        regroup_cell(&mut hits, group_buckets);
        assert_eq!(hits, again);
    }
    assert_eq!(q_off_bucket(entry_q_off(2)), 0);
    assert_eq!(q_off_bucket(entry_q_off(2049)), 0);
    assert_eq!(q_off_bucket(entry_q_off(2050)), 1);
    assert_eq!(entry_group(2 + 64 * 2048, 64), 1);
    assert_eq!(entry_group(1 + 64 * 2048, 64), 0);
}
