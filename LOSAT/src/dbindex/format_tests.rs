use super::*;

#[test]
fn a_written_header_reads_back_and_gives_ncbis_minimum_word_size() {
    for (legacy, minimum) in [(true, 12), (false, 14)] {
        let header = Header::of(legacy, 12, 3, 20, [7, 0, 19, 0]);
        let mut bytes = Vec::new();
        header.write(&mut bytes).unwrap();
        assert_eq!(bytes.len(), HEADER_SIZE);
        assert_eq!(bytes[0], if legacy { 5 } else { 6 });
        assert_eq!(Header::read(&bytes), Some(header));
        assert_eq!(header.min_word_size(), minimum);
    }
    // Past the end of a short file the mapping reads zeros; an Nmer over 15 fails.
    let short = [5u8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0];
    assert_eq!(Header::read(&short).unwrap().min_word_size(), 7);
    assert_eq!(Header::read(&[]).unwrap().min_word_size(), -1);
    let mut wide = short;
    wide[16] = 16;
    assert_eq!(Header::read(&wide), None);
}

#[test]
fn super_headers_and_volume_counts_read_as_ncbi() {
    let dir = crate::utils::test_dir::TempDir::new(&format!(
        "losat-dbindex-format-{}",
        crate::utils::test_dir::process_id()
    ));
    let shd = dir.join("v.shd");
    let header = SuperHeader {
        num_seq: 6,
        num_vol: 2,
    };
    assert!(!header.save(shd.to_str().unwrap().as_bytes()));
    assert_eq!(SuperHeader::load(&shd), Some(header));
    let mut bytes = std::fs::read(&shd).unwrap();
    bytes[4] = 2;
    std::fs::write(&shd, &bytes).unwrap();
    assert_eq!(SuperHeader::load(&shd), None);
    std::fs::write(&shd, &bytes[..8]).unwrap();
    assert_eq!(SuperHeader::load(&shd), None);
    assert_eq!(SuperHeader::load(&dir), None);
    let idx = dir.join("v.00.idx");
    let mut volume = Vec::new();
    Header::of(true, 8, 5, 28, [3, 0, 9, 0])
        .write(&mut volume)
        .unwrap();
    for (length, oids) in [
        (44, VolumeOids::Count(6)),
        (40, VolumeOids::Count(6)),
        (39, VolumeOids::Undefined),
        (36, VolumeOids::Undefined),
        (35, VolumeOids::ReadFailure("at end oid")),
        (28, VolumeOids::ReadFailure("at end oid")),
        (27, VolumeOids::ReadFailure("at start oid")),
        (0, VolumeOids::ReadFailure("at start oid")),
    ] {
        std::fs::write(&idx, &volume[..length]).unwrap();
        assert_eq!(volume_oids(&idx), oids, "{length}");
    }
    assert_eq!(volume_oids(&dir.join("none")), VolumeOids::Undefined);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn super_header_save_reports_closed_pipe_and_ignores_other_write_errors() {
    struct FailWriter(std::io::ErrorKind);
    impl std::io::Write for FailWriter {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::from(self.0))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::from(self.0))
        }
    }

    let header = SuperHeader {
        num_seq: 6,
        num_vol: 2,
    };
    let mut bytes = Vec::new();
    assert!(!header.save_to(&mut bytes));
    assert_eq!(bytes.len(), SUPER_HEADER_SIZE as usize);
    assert_eq!(words(&bytes), [0, 1, 6, 2]);

    let mut closed = FailWriter(std::io::ErrorKind::BrokenPipe);
    assert_eq!(header.save_to(&mut closed), cfg!(not(windows)));

    let mut full = FailWriter(std::io::ErrorKind::WriteZero);
    assert!(!header.save_to(&mut full));
}
