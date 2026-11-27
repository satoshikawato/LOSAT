use super::*;
use crate::seqid::defline_set::BlastDefline;
use crate::seqid::ParsedSeqId;

fn line(ids: Vec<ParsedSeqId>, taxid: Option<i32>) -> BlastDefline {
    BlastDefline {
        ids,
        taxid,
        ..BlastDefline::default()
    }
}

fn set(lines: Vec<BlastDefline>) -> DeflineSet {
    DeflineSet::from(lines)
}

/// A positive taxid of any line makes the Bioseq's source; else the first line's taxid decides
/// (probes ca4/ca5: the second line's 9606 passes; taxids_p5 OID 6: none fails).
#[test]
fn a_taxid_is_any_positive_one_else_the_first_lines() {
    let gi = |n| ParsedSeqId::Gi(n);
    assert!(has_taxid(
        0,
        &set(vec![line(vec![gi(1)], None), line(vec![gi(2)], Some(9606))])
    )
    .unwrap_or(false));
    assert!(!has_taxid(0, &set(vec![line(vec![gi(1)], None)])).unwrap_or(true));
    assert!(!has_taxid(0, &set(vec![line(vec![gi(1)], Some(0))])).unwrap_or(true));
    assert!(has_taxid(0, &set(vec![line(vec![gi(1)], Some(-3))])).unwrap_or(false));
    // Without a positive taxid, the first line whose first id of the first id's type matches it
    // and that has a taxid decides (`GetTaxidForSeqid`).
    assert!(has_taxid(
        0,
        &set(vec![line(vec![gi(1)], None), line(vec![gi(1)], Some(-3))])
    )
    .unwrap_or(false));
    assert!(!has_taxid(
        0,
        &set(vec![line(vec![gi(1)], None), line(vec![gi(2)], Some(-3))])
    )
    .unwrap_or(true));
    assert!(!has_taxid(
        0,
        &set(vec![
            line(vec![gi(1)], None),
            line(vec![gi(1)], Some(0)),
            line(vec![gi(1)], Some(-3))
        ])
    )
    .unwrap_or(true));
    // Only the line's first id of that type is compared.
    assert!(!has_taxid(
        0,
        &set(vec![
            line(vec![gi(1)], None),
            line(vec![gi(2), gi(1)], Some(-3))
        ])
    )
    .unwrap_or(true));
    // A first line without ids: NCBI reads a null reference.
    assert!(has_taxid(0, &set(vec![line(Vec::new(), None)])).is_err());
}

/// Approved exception `Q-D2b-OV=A` (`PD-LOSAT-NCBI-DEFECTS`, Owner 2026-11-25 17:54 JST):
/// `-ends` near `i32::MIN` widens `count - n` to `i64` in `Volume::ends` instead of overflowing
/// signed 32-bit `int`. Drives the real `Tests::run` / `Volume::ends` / `check_or_find` /
/// `test_oid` path on the `p5` fixture database (9 OIDs).
#[test]
fn widened_ends_near_int_min_jumps_past_volume() {
    use std::path::Path;

    use crate::blastdbcheck::log::{Log, MINUTIAE};
    use crate::blastdbcheck::TOOL;
    use crate::blastinput::ncbi_environment::DiagTarget;
    use crate::ncbi_app::files::standard_output_file;
    use crate::ncbi_app::Diag;
    use crate::report::run_files::{OpenFiles, RunFiles};
    use crate::seqdb::{SearchPath, SeqDb, SeqType};

    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/blastdb_aliastool/in");
    let path = SearchPath::new(&dir, None, None);
    let db = SeqDb::open_in(b"p5", SeqType::Protein, &path).unwrap();
    assert_eq!(db.num_oids(), 9);

    for (ends, expected_oid) in [
        (i32::MIN, 2_147_483_657i64),
        (-2_147_483_638i32, 2_147_483_647),
    ] {
        let mut tests = Tests::new(
            vec![Test::Ends(ends)],
            Flags {
                isam: false,
                legacy: false,
                taxids: false,
            },
            path.clone(),
        );
        let mut out = Vec::new();
        let mut err = Vec::new();
        let mut seen = Seen::default();
        let target = DiagTarget::Stderr;
        let faults = {
            let files = OpenFiles::new(&RunFiles::Disk);
            let mut log = Log::new(standard_output_file(&mut out, &files), MINUTIAE);
            let mut diag = Diag::new(&target, &mut err, &TOOL);
            match tests.run(&mut log, &db, b"p5", &mut seen, &mut diag) {
                Ok(faults) => faults,
                Err(Stop::Closed) => panic!("ends={ends}: unexpected closed stream"),
                Err(Stop::Refused(reason)) => panic!("ends={ends}: unexpected refusal: {reason}"),
            }
        };
        assert_eq!(faults, 1, "ends={ends}");
        assert!(
            err.is_empty(),
            "ends={ends}: {}",
            String::from_utf8_lossy(&err)
        );
        let text = String::from_utf8(out).unwrap();
        let expected = format!("Failed during headers, oid {expected_oid}: OID not found");
        assert!(text.contains(&expected), "ends={ends}:\n{text}");
    }
}
