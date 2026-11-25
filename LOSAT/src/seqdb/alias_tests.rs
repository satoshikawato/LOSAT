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

#[test]
fn v5_compare_seq_id_and_set_ops_match_ncbi_semantics() {
    use super::{compare_seq_id, merge_sort_by_result, set_difference_by, set_intersection_by};

    let cmp = |a: &[u8], b: &[u8]| compare_seq_id(a, b, b"test.txt");

    // Versioned vs unversioned: equal under s_CompareSeqId (neither is less than the other).
    assert!(!cmp(b"NP_000103.1", b"NP_000103").unwrap());
    assert!(!cmp(b"NP_000103", b"NP_000103.1").unwrap());
    // Distinct versions order lexicographically.
    assert!(cmp(b"NP_000103.1", b"NP_000103.2").unwrap());
    assert!(!cmp(b"NP_000103.2", b"NP_000103.1").unwrap());
    // Local ids compare case-insensitively under CSeq_id::Compare.
    assert!(!cmp(b"loc7", b"LOC7").unwrap());
    assert!(!cmp(b"LOC7", b"loc7").unwrap());
    // Unparseable FASTA form raises NCBI's CSeqIdException error.
    assert!(cmp(b"sp||", b"NP_000103.1").is_err());

    // Stable merge sort with compare_seq_id, then set intersection and difference.
    // Because compare_seq_id falls back to raw byte order for non-equivalent ids, "LOC7" sorts
    // before "NP_000103" while "loc7" sorts after "P10102.1" (matching NCBI's std::set_intersection
    // quirk in db.il.bp.v5.sq.slloc.6).
    let mut a = vec![
        b"NP_000118.1".to_vec(),
        b"NP_000103".to_vec(),
        b"LOC7".to_vec(),
    ];
    let mut b = vec![
        b"loc7".to_vec(),
        b"NP_000103.1".to_vec(),
        b"P10102.1".to_vec(),
    ];
    let mut less = |x: &[u8], y: &[u8]| compare_seq_id(x, y, b"test.txt");
    merge_sort_by_result(&mut a, &mut less).unwrap();
    merge_sort_by_result(&mut b, &mut less).unwrap();

    let inter = set_intersection_by(&a, &b, &mut less).unwrap();
    assert_eq!(inter, vec![b"NP_000103".to_vec()]);
    // Single-element case-different local ids compare equal directly.
    assert_eq!(
        set_intersection_by(&[b"LOC7".to_vec()], &[b"loc7".to_vec()], &mut less).unwrap(),
        vec![b"LOC7".to_vec()]
    );

    let diff = set_difference_by(&a, &b, &mut less).unwrap();
    assert_eq!(diff, vec![b"LOC7".to_vec(), b"NP_000118.1".to_vec()]);
}

#[test]
fn v5_alias_seqidlist_and_taxidlist_filter_oids_and_lines() {
    let idlist = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/idlist");
    let open_id = |name: &str, seqtype: SeqType| {
        SeqDb::open(idlist.join(name).to_str().unwrap(), seqtype).unwrap()
    };

    // v5 SEQIDLIST selects the same OIDs and keeps only matching deflines.
    let db_sq = open_id("g_sq", SeqType::Protein);
    assert_eq!(lines(&db_sq), [(0, 1), (1, 1), (6, 1), (13, 1)]);
    assert!(db_sq.drops_lineless_hits());

    // Binary .bsl SEQIDLIST selects the exact same OIDs and lines as text.
    let db_sqb = open_id("g_sqb", SeqType::Protein);
    assert_eq!(lines(&db_sqb), [(0, 1), (1, 1), (6, 1), (13, 1)]);

    // v5 TAXIDLIST selects OIDs with taxid 9606 and keeps only taxid-9606 deflines.
    let db_tx = open_id("g_tx", SeqType::Protein);
    let tx_oids: Vec<u32> = (0..db_tx.num_seqs()).map(|i| db_tx.oid(i).get()).collect();
    assert_eq!(tx_oids, [0, 1, 3, 5, 6, 7, 10, 13]);
    assert_eq!(db_tx.defline_set(0).unwrap().deflines.len(), 1);
    assert!(db_tx.drops_lineless_hits());

    // Combined SEQIDLIST + TAXIDLIST in one node ORs OIDs into Pass 2's f_bits and keeps lines
    // matching either list (OID 0 keeps both its P10102.1 and taxid-9606 lines).
    let db_sq_tx = open_id("g_sq_tx", SeqType::Protein);
    assert_eq!(
        lines(&db_sq_tx),
        [
            (0, 2),
            (1, 1),
            (3, 1),
            (5, 1),
            (6, 1),
            (7, 1),
            (10, 1),
            (13, 1)
        ]
    );

    // Version 5 TAXIDLIST over mismatched lookup/header taxids (lx): selects OIDs 1, 2, 3,
    // erases every line of OIDs 1 and 2, and marks them as dropped hits.
    let db_lx = open_id("lx_tx", SeqType::Protein);
    assert_eq!(lines(&db_lx), [(1, 0), (2, 0), (3, 1)]);
    assert!(!db_lx.keeps_line(0).unwrap());
    assert!(!db_lx.keeps_line(1).unwrap());
    assert!(db_lx.keeps_line(2).unwrap());

    // Nucleotide v5 SEQIDLIST (LMDB acc2oid is case-sensitive, so lowercase ab000207.1 misses)
    // and TAXIDLIST.
    let db_nsq = open_id("n_sq", SeqType::Nucleotide);
    assert_eq!(lines(&db_nsq), [(1, 1)]);
    let db_ntx = open_id("n_tx", SeqType::Nucleotide);
    assert_eq!(lines(&db_ntx), [(0, 1), (2, 1), (3, 1)]);
}
