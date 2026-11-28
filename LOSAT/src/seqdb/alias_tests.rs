//! Tests of the alias key `OID_MASK_TYPE` (line RG, series D, D1b-2) over the databases of
//! `LOSAT/tests/db_search/oidmask` (`make_oidmask.sh` there; NCBI BLAST+ 2.17.0's makeblastdb and
//! blastdb_aliastool): which definition lines are left of each sequence, per volume, by the last
//! OID list node above it.

use std::path::PathBuf;

use super::super::{DbUse, SeqDb, SeqDbError, SeqType};

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

fn memb_fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/memb_bit/db")
}

fn open_memb(name: &str, seqtype: SeqType) -> Result<SeqDb, SeqDbError> {
    SeqDb::open(memb_fixtures().join(name).to_str().unwrap(), seqtype)
}

const MB1_LINES: [(u32, usize); 6] = [(0, 1), (1, 2), (2, 0), (3, 0), (4, 0), (5, 0)];
const MEMB_CONFLICT: Option<&[u8]> = Some(b"MEMB_BIT error: conflicting bit found.");

#[test]
fn memb_bit_selects_defline_memberships_and_tracks_filters_header_sets() {
    let db0 = open_memb("mb_0", SeqType::Protein).unwrap();
    assert!(!db0.drops_lineless_hits());
    assert!(!db0.filters_header_sets(db0.location_of(0).unwrap()));
    assert_eq!(
        lines(&db0),
        [(0, 2), (1, 2), (2, 1), (3, 1), (4, 2), (5, 1)]
    );

    let db1 = open_memb("mb_1", SeqType::Protein).unwrap();
    assert!(db1.drops_lineless_hits());
    assert!(db1.filters_header_sets(db1.location_of(0).unwrap()));
    assert_eq!(lines(&db1), MB1_LINES);
    assert_eq!(
        (0..6)
            .map(|oid| db1.erases_a_line_of_oid(oid).unwrap())
            .collect::<Vec<_>>(),
        [true, false, true, true, true, true]
    );

    // OID 1 dl1 carries memberships {1, 2}: only front() (1) is tested, so bit 2 erases it.
    let db2 = open_memb("mb_2", SeqType::Protein).unwrap();
    assert_eq!(
        lines(&db2),
        [(0, 1), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)]
    );

    // Bits 31 and 32 (INT_MIN = 0x8000_0000) work on OID 4.
    let db31 = open_memb("mb_31", SeqType::Protein).unwrap();
    let db32 = open_memb("mb_32", SeqType::Protein).unwrap();
    assert_eq!(db31.defline_set(4).unwrap().deflines.len(), 1);
    assert_eq!(db32.defline_set(4).unwrap().deflines.len(), 1);
    assert_eq!(
        db32.defline_set(4).unwrap().deflines[0].title.as_deref(),
        Some(b"seq4_dl0 bit32".as_slice())
    );
    assert_eq!(
        db31.defline_set(4).unwrap().deflines[0].title.as_deref(),
        Some(b"seq4_dl1 bit31".as_slice())
    );
}

#[test]
fn memb_bit_conflicts_and_parse_order_match_ncbi() {
    // Parse order in ComputeMasks: SEQIDLIST before MEMB_BIT before TAXIDLIST.
    let err_si = open_memb("ord_si_mb", SeqType::Protein).unwrap_err();
    assert!(
        err_si.to_string().contains("has multiple SEQID lists"),
        "{err_si}"
    );
    let err_tx = open_memb("ord_mb_tx", SeqType::Protein).unwrap_err();
    assert!(
        err_tx
            .to_string()
            .contains("Cannot convert string 'abc' to Uint8"),
        "{err_tx}"
    );

    // Deferred conflict without NeedTotalsScan: open succeeds, build_oid_list fails.
    let db = open_memb("conflict_plain", SeqType::Protein).unwrap();
    let err = db.build_oid_list().unwrap_err();
    assert_eq!(err.ncbi_message().as_deref(), MEMB_CONFLICT);
    let site = super::super::throw_site(b"MEMB_BIT error: conflicting bit found.").unwrap();
    assert_eq!(
        (site.file, site.line, site.function),
        (
            "include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp",
            881,
            None
        )
    );

    // Parent 0 over child 1 conflicts; parent 1 over child 0 succeeds.
    let db_0_1 = open_memb("nest_0_1", SeqType::Protein).unwrap();
    assert_eq!(
        db_0_1
            .build_oid_list()
            .unwrap_err()
            .ncbi_message()
            .as_deref(),
        MEMB_CONFLICT
    );
    let db_1_0 = open_memb("nest_1_0", SeqType::Protein).unwrap();
    assert!(db_1_0.build_oid_list().is_ok());
    assert_eq!(lines(&db_1_0), MB1_LINES);

    // Conflict with NeedTotalsScan (GILIST without NSEQ+LENGTH): fails at open.
    let err_gi = open_memb("conflict_gi", SeqType::Protein).unwrap_err();
    assert_eq!(err_gi.ncbi_message().as_deref(), MEMB_CONFLICT);
}

#[test]
fn memb_bit_past_32_refuses_memberships_and_masklist_builds_oid_list() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-memb33-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let alias = dir.join("mb33.pal");
    std::fs::write(
        &alias,
        format!(
            "DBLIST {}\nMEMB_BIT 33\n",
            memb_fixtures().join("pm5").display()
        ),
    )
    .unwrap();
    let db = SeqDb::open(dir.join("mb33").to_str().unwrap(), SeqType::Protein).unwrap();
    // OID 2 has no memberships: deterministically erased without UB.
    assert_eq!(db.defline_set(2).unwrap().deflines.len(), 0);
    // OID 0 carries memberships: refused explicitly.
    let err = db.defline_set(0).unwrap_err();
    assert!(
        matches!(&err, SeqDbError::Unsupported { message, .. } if message.contains("past the 32 bits")),
        "{err}"
    );
    std::fs::remove_dir_all(&dir).unwrap();

    // MASKLIST + conflicting MEMB_BIT under NSEQ+LENGTH: read_mask_block reports OID-list error.
    let nog_cc = memb_fixtures().join("nogless").join("cc");
    let db_cc = SeqDb::open(nog_cc.to_str().unwrap(), SeqType::Nucleotide).unwrap();
    let range = db_cc.next_block(0, u64::MAX).unwrap();
    let mut block = super::super::MaskBlock::default();
    db_cc.read_mask_block(range, 0, &mut block);
    let err = block.ranges(0).unwrap_err();
    assert_eq!(err.error().ncbi_message().as_deref(), MEMB_CONFLICT);
}

#[test]
fn database_lookup_files_collects_alias_index_group_and_existing_dblist() {
    use super::database_lookup_files;
    // Plain nucleotide and protein names yield group alias + alias + index in lookup order.
    assert_eq!(
        database_lookup_files(b"x", SeqType::Nucleotide, None, DbUse::search("blastn")),
        vec![b"index.alx".to_vec(), b"x.nal".to_vec(), b"x.nin".to_vec()]
    );
    assert_eq!(
        database_lookup_files(b"p5", SeqType::Protein, None, DbUse::search("blastn")),
        vec![
            b"index.alx".to_vec(),
            b"p5.pal".to_vec(),
            b"p5.pin".to_vec()
        ]
    );
    // Multiple quoted names with spaces: stops at the first missing name (matching NCBI's
    // `resolve_in` short-circuit), and continues to the second name when the first exists.
    assert_eq!(
        database_lookup_files(
            b"\"sub dir/x1\" x2",
            SeqType::Nucleotide,
            None,
            DbUse::search("blastn")
        ),
        vec![
            b"sub dir/index.alx".to_vec(),
            b"sub dir/x1.nal".to_vec(),
            b"sub dir/x1.nin".to_vec(),
        ]
    );
    let qdir = crate::utils::test_dir::TempDir::unique("losat-d19b-quoted");
    let sub = qdir.join("sub dir");
    std::fs::create_dir_all(&sub).unwrap();
    std::fs::write(sub.join("x1.nin"), b"x").unwrap();
    let quoted_arg = format!("\"{}\" x2", sub.join("x1").display());
    let quoted_got = database_lookup_files(
        quoted_arg.as_bytes(),
        SeqType::Nucleotide,
        None,
        DbUse::search("blastn"),
    );
    assert!(quoted_got.iter().any(|p| p == b"x2.nal"), "{quoted_got:?}");
    assert!(quoted_got.iter().any(|p| p == b"x2.nin"), "{quoted_got:?}");
    // Existing alias file expands DBLIST and stops on mutual recursion (before `compute_masks`
    // runs, so filter list paths like `SEQIDLIST` are not recorded on resolution failure).
    let dir = crate::utils::test_dir::TempDir::unique("losat-d19b-alias");
    std::fs::write(
        dir.join("a.nal"),
        b"TITLE a\nDBLIST v1 b\nSEQIDLIST ids.txt\n",
    )
    .unwrap();
    std::fs::write(dir.join("v1.nin"), b"x").unwrap();
    std::fs::write(dir.join("b.nal"), b"TITLE b\nDBLIST a\n").unwrap();
    let base = dir.join("a");
    let got = database_lookup_files(
        base.as_os_str().as_encoded_bytes(),
        SeqType::Nucleotide,
        None,
        DbUse::search("blastn"),
    );
    let has = |list: &[Vec<u8>], path: &std::path::Path| {
        let want = path.as_os_str().as_encoded_bytes();
        list.iter().any(|p| p.as_slice() == want)
    };
    assert!(has(&got, &dir.join("a.nal")));
    assert!(has(&got, &dir.join("b.nal")));
    assert!(has(&got, &dir.join("v1.nal")));
    assert!(has(&got, &dir.join("v1.nin")));
    assert!(!has(&got, &dir.join("ids.txt")));
    assert!(has(&got, &dir.join("index.alx")));

    // Valid alias records `compute_masks` paths: `SEQIDLIST`, `OIDLIST` primary + first-volume
    // fallback when primary is missing, and (on Unix) preserves backslashes literally just like
    // `compute_masks` (without `fix_delimiters`).
    let mask_dir = crate::utils::test_dir::TempDir::unique("losat-d19b-masks");
    let vol_dir = mask_dir.join("vols");
    std::fs::create_dir_all(&vol_dir).unwrap();
    std::fs::write(vol_dir.join("v1.nin"), b"x").unwrap();
    std::fs::write(
        mask_dir.join("m.nal"),
        b"TITLE m\nDBLIST v1\nOIDLIST mask.msk\nSEQIDLIST ids.txt\n",
    )
    .unwrap();
    let lookup_mask = |alias_name: &str| {
        database_lookup_files(
            mask_dir.join(alias_name).as_os_str().as_encoded_bytes(),
            SeqType::Nucleotide,
            Some(vol_dir.as_os_str().as_encoded_bytes()),
            DbUse::search("blastn"),
        )
    };
    let got_m = lookup_mask("m");
    assert!(has(&got_m, &mask_dir.join("mask.msk")), "{got_m:?}");
    assert!(has(&got_m, &vol_dir.join("mask.msk")), "{got_m:?}");
    assert!(has(&got_m, &mask_dir.join("ids.txt")), "{got_m:?}");
    // When primary OIDLIST already exists, `compute_masks` does not consult the volume fallback.
    std::fs::write(mask_dir.join("mask.msk"), b"x").unwrap();
    let got_m_present = lookup_mask("m");
    assert!(has(&got_m_present, &mask_dir.join("mask.msk")));
    assert!(
        !has(&got_m_present, &vol_dir.join("mask.msk")),
        "{got_m_present:?}"
    );
    #[cfg(unix)]
    {
        std::fs::write(
            mask_dir.join("bs.nal"),
            b"TITLE bs\nDBLIST v1\nOIDLIST sub\\mask.msk\nSEQIDLIST sub\\ids.txt\n",
        )
        .unwrap();
        let got_bs = lookup_mask("bs");
        assert!(has(&got_bs, &mask_dir.join(r"sub\mask.msk")), "{got_bs:?}");
        assert!(has(&got_bs, &vol_dir.join(r"sub\mask.msk")), "{got_bs:?}");
        assert!(has(&got_bs, &mask_dir.join(r"sub\ids.txt")), "{got_bs:?}");
        assert!(!has(&got_bs, &mask_dir.join("sub/mask.msk")), "{got_bs:?}");
        assert!(!has(&got_bs, &mask_dir.join("sub/ids.txt")), "{got_bs:?}");
    }

    // Finding 2 (a): local copy in cwd (`!skip_local`) is checked for DBLIST entries without a
    // directory when the alias sits in a BLASTDB directory.
    let case_b = crate::utils::test_dir::TempDir::unique("losat-d19b-local-copy");
    let b_dbs = case_b.join("dbs");
    std::fs::create_dir_all(&b_dbs).unwrap();
    std::fs::write(b_dbs.join("a.pal"), b"TITLE a\nDBLIST v1\n").unwrap();
    std::fs::write(b_dbs.join("v1.pin"), b"x").unwrap();
    let got_b = database_lookup_files(
        b"a",
        SeqType::Protein,
        Some(b_dbs.as_os_str().as_encoded_bytes()),
        DbUse::search("blastn"),
    );
    assert!(
        got_b.iter().any(|p| p == b"v1.pin"),
        "local copy v1.pin in cwd must be recorded: {got_b:?}"
    );

    // Finding 2 (b): DBLIST entry resolved via BLASTDB search path (and recursive alias inside
    // BLASTDB) records every checked `.pal`/`.pin` along the search path.
    let case_c = crate::utils::test_dir::TempDir::unique("losat-d19b-search-relookup");
    let c_dbs = case_c.join("dbs");
    std::fs::create_dir_all(&c_dbs).unwrap();
    std::fs::write(case_c.join("a.pal"), b"TITLE a\nDBLIST sub\n").unwrap();
    std::fs::write(c_dbs.join("sub.pal"), b"TITLE sub\nDBLIST v1\n").unwrap();
    std::fs::write(c_dbs.join("v1.pin"), b"x").unwrap();
    let got_c = database_lookup_files(
        case_c.join("a").as_os_str().as_encoded_bytes(),
        SeqType::Protein,
        Some(c_dbs.as_os_str().as_encoded_bytes()),
        DbUse::search("blastn"),
    );
    assert!(has(&got_c, &c_dbs.join("sub.pal")), "{got_c:?}");
    assert!(has(&got_c, &c_dbs.join("v1.pal")), "{got_c:?}");
    assert!(has(&got_c, &c_dbs.join("v1.pin")), "{got_c:?}");
}

// Line WQA D26: LOSAT Web's files of a database are those that the command line refuses
// (`SeqDb::read_paths`): besides the lookup candidates and the alias lists, each volume's
// component files (missing ones too) and the `MASKLIST` GI-mask files; the `OIDLIST` path beside
// the first volume only where the one beside the alias file is not a file.
#[test]
fn database_lookup_files_hold_the_volume_and_gi_mask_files() {
    use super::database_lookup_files;
    let dir = crate::utils::test_dir::TempDir::unique("losat-d26-lookup-files");
    let vols = dir.join("vols");
    std::fs::create_dir_all(&vols).unwrap();
    for ext in ["nin", "nhr", "nsq"] {
        std::fs::write(vols.join(format!("v.{ext}")), b"x").unwrap();
    }
    std::fs::write(dir.join("o.msk"), b"x").unwrap();
    std::fs::write(
        dir.join("m.nal"),
        format!(
            "TITLE m\nDBLIST v\nOIDLIST o.msk\nMASKLIST {}\n",
            dir.join("mk").display()
        ),
    )
    .unwrap();
    let got = database_lookup_files(
        dir.join("m").as_os_str().as_encoded_bytes(),
        SeqType::Nucleotide,
        Some(vols.as_os_str().as_encoded_bytes()),
        DbUse::TOOL,
    );
    let has = |path: &std::path::Path| {
        let want = path.as_os_str().as_encoded_bytes();
        got.iter().any(|p| p.as_slice() == want)
    };
    for file in [
        dir.join("m.nal"),
        dir.join("o.msk"),
        dir.join("index.alx"),
        vols.join("v.nin"),
        vols.join("v.nsq"),
        vols.join("v.nog"),
        vols.join("v.ndb"),
        dir.join("mk.gni"),
        dir.join("mk.gno"),
    ] {
        assert!(has(&file), "{} in {got:?}", file.display());
    }
    for file in [vols.join("o.msk"), dir.join("mk.gnd"), vols.join("v.pin")] {
        assert!(!has(&file), "{} in {got:?}", file.display());
    }
    // The files whose presence the resolution checks come apart from those of what it found
    // (`windowmasker -out` counts only the former).
    let split = super::database_files(
        dir.join("m").as_os_str().as_encoded_bytes(),
        SeqType::Nucleotide,
        Some(vols.as_os_str().as_encoded_bytes()),
        DbUse::TOOL,
    );
    let within = |files: &[Vec<u8>], path: std::path::PathBuf| {
        files
            .iter()
            .any(|p| p.as_slice() == path.as_os_str().as_encoded_bytes())
    };
    for file in [dir.join("m.nal"), dir.join("o.msk"), dir.join("index.alx")] {
        assert!(within(&split.looked_up, file.clone()), "{}", file.display());
        assert!(!within(&split.found, file.clone()), "{}", file.display());
    }
    for file in [vols.join("v.nog"), vols.join("v.ndb"), dir.join("mk.gni")] {
        assert!(within(&split.found, file.clone()), "{}", file.display());
        assert!(
            !within(&split.looked_up, file.clone()),
            "{}",
            file.display()
        );
    }
    assert_eq!(split.looked_up.len() + split.found.len(), got.len());
}

// Line WQA D26 round 1: the files of a database that a run can open follow its program and its
// options (`DbUse`): a search of the sequences opens the volume's index, header, sequence and
// string-index files and the LMDB file; the id-lookup files only with an id list or a taxonomy
// filter, the column files and GI-mask files only with a mask option, the RPS files only for the
// RPS searches; a name that none of these can open is not a file the run reads.
#[test]
fn database_files_follow_what_the_run_asks_for() {
    use super::database_lookup_files;
    let dir = crate::utils::test_dir::TempDir::unique("losat-d26-use");
    for ext in ["nin", "nhr", "nsq", "naa"] {
        std::fs::write(dir.join(format!("v.{ext}")), b"x").unwrap();
    }
    for ext in ["pin", "phr", "psq"] {
        std::fs::write(dir.join(format!("w.{ext}")), b"x").unwrap();
    }
    std::fs::write(
        dir.join("m.nal"),
        format!("TITLE m\nDBLIST v\nMASKLIST {}\n", dir.join("mk").display()),
    )
    .unwrap();
    let files = |name: &[u8], seqtype: SeqType, db_use: DbUse| -> Vec<String> {
        let base = dir.as_os_str().as_encoded_bytes();
        let name = [base, b"/", name].concat();
        database_lookup_files(&name, seqtype, Some(base), db_use)
            .iter()
            .map(|file| {
                let file = String::from_utf8_lossy(file).into_owned();
                file.strip_prefix(&format!("{}/", dir.display()))
                    .map_or(file.clone(), str::to_string)
            })
            .collect()
    };
    let with = |option: &str| {
        let mut db_use = DbUse::search("blastn");
        db_use.note_option(option);
        db_use
    };
    let has = |files: &[String], name: &str| files.iter().any(|file| file == name);
    // The plain search.
    let plain = files(b"m", SeqType::Nucleotide, DbUse::search("blastn"));
    for name in [
        "m.nal", "v.nin", "v.nhr", "v.nsq", "v.nsd", "v.nsi", "v.ndb",
    ] {
        assert!(has(&plain, name), "{name} in {plain:?}");
    }
    for name in [
        "mk.gni", "mk.gno", "v.nog", "v.naa", "v.nba", "v.nnd", "v.nni", "v.nos", "v.nto", "v.njs",
        "v.nhd", "v.pin", "v.rps", "m.njs",
    ] {
        assert!(!has(&plain, name), "{name} in {plain:?}");
    }
    // A mask option: the GI masks, the GI index, and each volume's columns (the data files only
    // where the column's index file exists).
    for option in ["db_soft_mask", "db_hard_mask"] {
        let masked = files(b"m", SeqType::Nucleotide, with(option));
        for name in [
            "mk.gni", "mk.gno", "v.nog", "v.naa", "v.nab", "v.nac", "v.nba", "v.nza",
        ] {
            assert!(has(&masked, name), "{option}: {name} in {masked:?}");
        }
        for name in ["v.nbb", "v.nbc", "v.nnd", "v.nos"] {
            assert!(!has(&masked, name), "{option}: {name} in {masked:?}");
        }
    }
    // An imported search strategy may hold any of the options, so it asks for all of them.
    let imported = files(b"m", SeqType::Nucleotide, with("import_search_strategy"));
    for name in [
        "mk.gni", "v.nog", "v.naa", "v.nnd", "v.nos", "v.nto", "v.nti",
    ] {
        assert!(has(&imported, name), "{name} in {imported:?}");
    }
    // An id list or a taxonomy filter: the lookup files of its kind, not the masks nor the
    // lookup files of another kind.
    for (option, asked, others) in [
        (
            "gilist",
            &["v.nni", "v.nnd"][..],
            &["v.nos", "v.nto", "v.nti", "v.npi"][..],
        ),
        (
            "negative_gilist",
            &["v.nni", "v.nnd"],
            &["v.nos", "v.nto", "v.nti"],
        ),
        ("ipglist", &["v.npi", "v.npd"], &["v.nni", "v.nos", "v.nto"]),
        (
            "seqidlist",
            &["v.nos", "v.nti", "v.ntd"],
            &["v.nni", "v.nnd", "v.nto", "v.ntf"],
        ),
        (
            "negative_seqidlist",
            &["v.nos", "v.nti", "v.ntd"],
            &["v.nni", "v.nto"],
        ),
        (
            "taxids",
            &["v.not", "v.ntf", "v.nto"],
            &["v.nni", "v.nos", "v.nti"],
        ),
        (
            "negative_taxids",
            &["v.not", "v.ntf", "v.nto"],
            &["v.nni", "v.nos"],
        ),
        (
            "taxidlist",
            &["v.not", "v.ntf", "v.nto"],
            &["v.nni", "v.nos"],
        ),
        (
            "negative_taxidlist",
            &["v.not", "v.ntf", "v.nto"],
            &["v.nni", "v.nos"],
        ),
    ] {
        let filtered = files(b"m", SeqType::Nucleotide, with(option));
        for name in asked {
            assert!(has(&filtered, name), "{option}: {name} in {filtered:?}");
        }
        for name in others
            .iter()
            .chain(&["mk.gni", "v.nog", "v.naa", "v.nhd", "v.njs"])
        {
            assert!(!has(&filtered, name), "{option}: {name} in {filtered:?}");
        }
    }
    // The RPS files belong to the RPS searches only.
    let blastp = files(b"w", SeqType::Protein, DbUse::search("blastp"));
    let rps = files(b"w", SeqType::Protein, DbUse::search("rpsblast"));
    for name in ["w.pin", "w.phr", "w.psq", "w.pdb"] {
        assert!(has(&blastp, name) && has(&rps, name), "{name}");
    }
    for name in ["w.rps", "w.loo", "w.aux", "w.freq"] {
        assert!(!has(&blastp, name), "{name} in {blastp:?}");
        assert!(has(&rps, name), "{name} in {rps:?}");
    }
    // A tool that is not a search may open all of them, but not the RPS files.
    let tool = files(b"m", SeqType::Nucleotide, DbUse::TOOL);
    for name in [
        "v.nog", "v.naa", "v.nnd", "v.nos", "v.nhd", "v.njs", "mk.gni", "m.njs",
    ] {
        assert!(has(&tool, name), "{name} in {tool:?}");
    }
}

#[test]
fn v5_compare_seq_id_and_set_ops_match_ncbi_semantics() {
    use super::{compare_seq_id, set_difference_by, set_intersection_by, validate_compare_seqids};

    let cmp = |a: &[u8], b: &[u8]| {
        compare_seq_id(
            a,
            Some(b"test.txt"),
            "an alias file's SEQIDLIST",
            b,
            None,
            "the search's -seqidlist",
        )
    };

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

    // Lists without differently spelled equivalent ids pass the strict-weak-ordering guard and
    // sort uniquely by bytewise `<` (matching NCBI's std::sort + set_intersection/set_difference).
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
    a.sort_unstable();
    b.sort_unstable();
    validate_compare_seqids(a.as_slice(), Some(b"test.txt"), "an alias file's SEQIDLIST").unwrap();
    validate_compare_seqids(b.as_slice(), None, "the search's -seqidlist").unwrap();

    let mut less_lr = |x: &[u8], y: &[u8]| cmp(x, y);
    let mut less_rl = |x: &[u8], y: &[u8]| cmp(x, y);
    let inter =
        set_intersection_by(a.as_slice(), b.as_slice(), &mut less_lr, &mut less_rl).unwrap();
    assert_eq!(inter, vec![b"NP_000103".as_slice()]);
    // Single-element case-different local ids compare equal directly.
    let single_lhs = [b"LOC7".to_vec()];
    let single_rhs = [b"loc7".to_vec()];
    assert_eq!(
        set_intersection_by(
            single_lhs.as_slice(),
            single_rhs.as_slice(),
            &mut less_lr,
            &mut less_rl
        )
        .unwrap(),
        vec![b"LOC7".as_slice()]
    );

    let diff = set_difference_by(a.as_slice(), b.as_slice(), &mut less_lr, &mut less_rl).unwrap();
    assert_eq!(diff, vec![b"LOC7".as_slice(), b"NP_000118.1".as_slice()]);

    // Multiset semantics of std::set_difference: duplicate entries are preserved when not deduped.
    let dup_lhs = [b"NP_000103.1".to_vec(), b"NP_000103.1".to_vec()];
    let dup_rhs = [b"NP_000103.1".to_vec()];
    let dup_diff = set_difference_by(
        dup_lhs.as_slice(),
        dup_rhs.as_slice(),
        &mut less_lr,
        &mut less_rl,
    )
    .unwrap();
    assert_eq!(dup_diff, vec![b"NP_000103.1".as_slice()]);
}

#[test]
fn v5_rejects_non_strict_weak_equivalent_seqid_pairs() {
    use super::validate_compare_seqids;

    let check = |ids: &[&[u8]]| {
        let mut owned: Vec<Vec<u8>> = ids.iter().map(|s| s.to_vec()).collect();
        owned.sort_unstable();
        validate_compare_seqids(
            owned.as_slice(),
            Some(b"test.txt"),
            "an alias file's SEQIDLIST",
        )
    };

    // Unversioned + versioned pairs (all three permutations from the review probe), case-variant
    // local/accession pairs, cross-database accession pairs, and chainless/chained PDB pairs must
    // be explicitly refused.
    for pair in [
        &[
            b"NP_000103".as_slice(),
            b"NP_000103.1".as_slice(),
            b"NP_000103.2".as_slice(),
        ][..],
        &[b"NP_000103.2", b"NP_000103", b"NP_000103.1"],
        &[b"NP_000103.1", b"NP_000103.2", b"NP_000103"],
        &[b"LOC7", b"loc7"],
        &[b"NP_000103.1", b"np_000103.1"],
        &[b"gb|AB000001.1|", b"emb|AB000001.1|"],
        &[b"pdb|1ABC|A", b"pdb|1ABC|"],
        &[b"123", b"gi|123"],
        &[b"loc7", b"lcl|loc7"],
    ] {
        let err = check(pair).unwrap_err();
        let text = err.to_error().to_string();
        assert!(
            text.contains("holds two differently spelled Seq-ids"),
            "unexpected error for {pair:?}: {text}"
        );
    }

    // Distinct versions without the unversioned form, unsorted non-equivalent ids, bare GI vs
    // local integer (`123` parses as GI under `fParse_AnyRaw`), single-entry lists, and identical
    // bytewise duplicates (even unparseable ones, where s_CompareSeqId short-circuits id1 == id2)
    // are allowed.
    for ok in [
        &[][..],
        &[b"NP_000103".as_slice()],
        &[b"NP_000118.1", b"NP_000103.1"],
        &[b"NP_000103.2", b"NP_000103.1", b"LOC7"],
        &[b"123", b"lcl|123"],
        &[b"NP_000103.1", b"NP_000103.1"],
        &[b"sp||", b"sp||"],
    ] {
        assert!(check(ok).is_ok(), "unexpected refusal for {ok:?}");
    }

    // One common exception text across all unparseable entries raises NCBI's CSeqIdException
    // (including two distinct NUL-prefixed tokens whose C-string `what()` is identical).
    for unparseable in [
        &[b"sp||".as_slice(), b"NP_000103.1".as_slice()][..],
        &[b"NP_000103", b"NP_000103.1", b"sp||"],
        &[b"\x00aaa", b"\x00bbb", b"NP_000103.1"],
    ] {
        let err = check(unparseable).unwrap_err().to_error();
        assert!(err.ncbi_exception().is_some(), "{unparseable:?}: {err}");
    }

    // Two or more distinct exception texts depend on the C++ library's sort order and are refused
    // identically across all input permutations.
    for multi_exc in [
        &[
            b"sp||".as_slice(),
            b"tr||".as_slice(),
            b"NP_000103.1".as_slice(),
        ][..],
        &[b"tr||", b"sp||", b"NP_000103.1"],
        &[b"NP_000103.1", b"sp||", b"tr||"],
        &[b"NP_000103.1", b"tr||", b"sp||"],
        &[b"\x00aaa", b"PM\x00bbb"],
    ] {
        let err = check(multi_exc).unwrap_err().to_error().to_string();
        assert!(
            err.contains("fail NCBI's `s_CompareSeqId` with different exceptions"),
            "{multi_exc:?}: {err}"
        );
    }

    // A Seq-id that posts a warning alongside an unparseable token is refused (step 1 before
    // step 3) regardless of whether the warning sorts before or after the unparseable token.
    // `AAAAA01P000001` (5+8 with protein flag `P` and division 0) posts accguide warning 11.
    for warn_and_bad in [
        [b"AAAAA01P000001".as_slice(), b"sp||".as_slice()],
        [b"\x00bad".as_slice(), b"AAAAA01P000001".as_slice()],
    ] {
        let err = check(&warn_and_bad).unwrap_err().to_error().to_string();
        assert!(
            err.contains("posts a Seq-id warning"),
            "{warn_and_bad:?}: {err}"
        );
    }

    // End-to-end SeqDb::restrict_seqids checks over the idlist fixtures:
    // 1. Without a search seqid list, an alias SEQIDLIST with equivalent spellings opens normally.
    let idlist = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/idlist");
    let open_id =
        |name: &str| SeqDb::open(idlist.join(name).to_str().unwrap(), SeqType::Protein).unwrap();
    let read_list = |name: &str, negative: bool| {
        let bytes = std::fs::read(idlist.join(name)).unwrap();
        let mut warns = Vec::new();
        let file = super::super::seqidlist::read_si_list(&bytes, &mut warns).unwrap();
        super::super::filter::SeqidList { negative, file }
    };
    for alias in ["g_sqamb0", "g_sqamb1", "g_sqamb2"] {
        let db = open_id(alias);
        assert_eq!(lines(&db), [(1, 1)], "{alias}");
        // 2. With -seqidlist or -negative_seqidlist, the ambiguous alias SEQIDLIST is refused.
        for (list_name, negative) in [
            ("u_np103v.txt", false),
            ("u_np103_1.txt", false),
            ("u_np103_2.txt", false),
            ("u_p10102.txt", true),
        ] {
            let mut db = open_id(alias);
            let mut warnings = Vec::new();
            let err = db
                .restrict_seqids(&read_list(list_name, negative), &mut warnings)
                .unwrap_err()
                .to_string();
            assert!(
                err.contains("an alias file's SEQIDLIST holds two differently spelled Seq-ids"),
                "{alias} + {list_name} (neg={negative}): {err}"
            );
        }
    }
    // 3. Clean alias SEQIDLIST + ambiguous -seqidlist is refused without naming the alias file
    // path, while ambiguous -negative_seqidlist succeeds (NCBI sorts negative_seqidlist with `<`).
    let mut db = open_id("g_sq");
    let mut warnings = Vec::new();
    let err = db
        .restrict_seqids(&read_list("u_np103_amb.txt", false), &mut warnings)
        .unwrap_err()
        .to_string();
    assert!(
        err.starts_with(
            "unsupported database: the search's -seqidlist holds two differently spelled Seq-ids"
        ),
        "{err}"
    );
    assert!(
        !err.contains("sq_v5.txt"),
        "must not name alias path: {err}"
    );
    let mut db = open_id("g_sq");
    let mut warnings = Vec::new();
    db.restrict_seqids(&read_list("u_np103_amb.txt", true), &mut warnings)
        .unwrap();
    assert_eq!(lines(&db), [(0, 1), (6, 1), (13, 1)]);
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
