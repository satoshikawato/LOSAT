//! Tests of the alias key `OID_MASK_TYPE` (line RG, series D, D1b-2) over the databases of
//! `LOSAT/tests/db_search/oidmask` (`make_oidmask.sh` there; NCBI BLAST+ 2.17.0's makeblastdb and
//! blastdb_aliastool): which definition lines are left of each sequence, per volume, by the last
//! OID list node above it.

use std::path::PathBuf;

use super::super::{SeqDb, SeqDbError, SeqType};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/oidmask")
}

fn open(name: &str, seqtype: SeqType) -> Result<SeqDb, SeqDbError> {
    SeqDb::open(fixtures().join(name).to_str().unwrap(), seqtype)
}

/// The OID and the number of definition lines left of each sequence the database includes.
fn lines(db: &SeqDb) -> Vec<(u32, usize)> {
    (0..db.num_seqs())
        .map(|at| (db.oid(at).get(), db.defline_set(at).unwrap().deflines.len()))
        .collect()
}

/// The protein records (`om.faa`): OIDs 1, 2 and 5 have predicted RefSeq ids alone (2 with a
/// GenBank id beside it), OID 3 two lines, the first predicted; `omp.pxm` leaves out 1 and 5.
const ALL_LINES: [(u32, usize); 8] = [
    (0, 1),
    (1, 1),
    (2, 1),
    (3, 2),
    (4, 1),
    (5, 1),
    (6, 1),
    (7, 1),
];
const MASKED_LINES: [(u32, usize); 6] = [(0, 1), (2, 1), (3, 2), (4, 1), (6, 1), (7, 1)];
const MASKED_ERASED: [(u32, usize); 6] = [(0, 1), (2, 0), (3, 1), (4, 1), (6, 1), (7, 1)];

#[test]
fn bit_zero_of_the_type_erases_the_lines_of_predicted_ids() {
    // NCBI's own alias (blastdb_aliastool -oid_masks 1), versions 5 and 4.
    for name in ["ompx", "omp4x"] {
        let db = open(name, SeqType::Protein).unwrap();
        assert_eq!(lines(&db), MASKED_ERASED, "{name}");
        assert!(db.drops_lineless_hits(), "{name}");
        let kept: Vec<bool> = (0..db.num_seqs())
            .map(|at| db.keeps_line(at).unwrap())
            .collect();
        assert_eq!(kept, [true, false, true, true, true, true], "{name}");
    }
    // Every OID kept: three are left without a line.
    let db = open("tall", SeqType::Protein).unwrap();
    assert_eq!(
        lines(&db),
        [
            (0, 1),
            (1, 0),
            (2, 0),
            (3, 1),
            (4, 1),
            (5, 0),
            (6, 1),
            (7, 1)
        ]
    );
    // The values NCBI's NStr::StringToUInt reads with bit 0 set, with the scan of the
    // constructor (no NSEQ and LENGTH) too.
    for name in ["t3", "tmax", "tplus", "tscan"] {
        assert_eq!(
            lines(&open(name, SeqType::Protein).unwrap()),
            MASKED_ERASED,
            "{name}"
        );
    }
    // Bit 0 clear, no key, or the key on a node without OIDLIST (not read, even as `abc`).
    for name in ["t0", "t2", "tnone"] {
        let db = open(name, SeqType::Protein).unwrap();
        assert_eq!(lines(&db), MASKED_LINES, "{name}");
        assert!(!db.drops_lineless_hits(), "{name}");
    }
    assert_eq!(lines(&open("tnol", SeqType::Protein).unwrap()), ALL_LINES);
}

#[test]
fn a_type_that_does_not_convert_stops_the_open() {
    for (name, text) in [
        ("tabc", "Cannot convert string 'abc' to Uint8"),
        (
            "tover",
            "Cannot convert string '4294967296' to unsigned int",
        ),
        ("tneg", "Cannot convert string '-1' to Uint8"),
        ("tempty", "Cannot convert string '' to Uint8"),
        ("tfirst", "Cannot convert string 'abc' to Uint8"),
    ] {
        let error = open(name, SeqType::Protein).err().unwrap().to_string();
        assert!(error.contains(text), "{name}: {error}");
    }
}

#[test]
fn the_last_oid_list_node_above_a_volume_decides() {
    // FIRST_OID 2 and LAST_OID 6 (the OIDs before it) with the list and the type on one node
    // (row db.om.bp.j1.6: NCBI reports OIDs 3 and 4).
    assert_eq!(
        lines(&open("j1", SeqType::Protein).unwrap()),
        [(2, 0), (3, 1), (4, 1)]
    );
    // A parent's OIDLIST without the key over ompx: the parent is evaluated last, type 0.
    assert_eq!(lines(&open("j2", SeqType::Protein).unwrap()), MASKED_LINES);
    // Two siblings over one volume: the later one's type, 0 (ompy) here, 1 (ompx) the other way.
    assert_eq!(lines(&open("j5", SeqType::Protein).unwrap()), ALL_LINES);
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-oidmask-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let reversed = dir.join("j5r.pal");
    let fixtures = fixtures();
    let alias = format!(
        "DBLIST {} {}\n",
        fixtures.join("ompy").display(),
        fixtures.join("ompx").display()
    );
    std::fs::write(&reversed, alias).unwrap();
    let db = SeqDb::open(dir.join("j5r").to_str().unwrap(), SeqType::Protein);
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(
        lines(&db.unwrap()),
        [
            (0, 1),
            (1, 0),
            (2, 0),
            (3, 1),
            (4, 1),
            (5, 0),
            (6, 1),
            (7, 1)
        ]
    );
}

#[test]
fn each_volume_has_the_type_of_its_own_list() {
    // omnx: an alias per volume (OIDLIST of its .nxm, OID_MASK_TYPE 1); ntall: one list and type
    // over both volumes.
    let db = open("omnx", SeqType::Nucleotide).unwrap();
    assert_eq!(lines(&db), MASKED_ERASED);
    let db = open("ntall", SeqType::Nucleotide).unwrap();
    assert_eq!(
        lines(&db),
        [
            (0, 1),
            (1, 0),
            (2, 0),
            (3, 1),
            (4, 1),
            (5, 0),
            (6, 1),
            (7, 1)
        ]
    );
    assert_eq!(lines(&open("omn", SeqType::Nucleotide).unwrap()), ALL_LINES);
}

/// `omb` holds the byte 0x80 in the titles of two predicted lines (`make_oidmask.sh`): NCBI
/// writes its error for each such byte at every read of the header, also in a line it then
/// erases (probe b1 of line RG D1b-2), which LOSAT refuses at the read; a sequence without such
/// a line is read.
#[test]
fn a_removed_line_with_bytes_ncbi_reads_as_hash_is_refused() {
    let db = open("ombx", SeqType::Protein).unwrap();
    let error = db.defline_set(0).unwrap_err().to_string();
    assert!(
        error.contains("removes, with bytes NCBI reads as '#'"),
        "{error}"
    );
    assert!(db.keeps_line(0).is_err());
    assert_eq!(db.defline_set(1).unwrap().deflines.len(), 1);
}

/// The reads of a volume that NCBI's RPS search of several volumes searches bare: the stored
/// header keeps the lines the type erases, and `erases_a_line_of_oid` tells which sequences lose
/// one (OIDs 1 and 5, which the list leaves out, too).
#[test]
fn the_stored_header_keeps_the_lines_the_type_erases() {
    use super::super::SequenceHeader;
    let db = open("ompx", SeqType::Protein).unwrap();
    let count = |header: SequenceHeader| match header {
        SequenceHeader::Parsed(set) => set.deflines.len(),
        SequenceHeader::Ordinal(deflines) => deflines.len(),
    };
    let erased: Vec<bool> = (0..8)
        .map(|oid| db.erases_a_line_of_oid(oid).unwrap())
        .collect();
    assert_eq!(erased, [false, true, true, true, false, true, false, false]);
    assert_eq!(count(db.read_header_of_oid(2).unwrap()), 0);
    assert_eq!(count(db.read_stored_header_of_oid(2).unwrap()), 1);
    assert_eq!(count(db.read_header_of_oid(3).unwrap()), 1);
    assert_eq!(count(db.read_stored_header_of_oid(3).unwrap()), 2);
    let db = open("ompy", SeqType::Protein).unwrap();
    assert!((0..8).all(|oid| !db.erases_a_line_of_oid(oid).unwrap()));
}
