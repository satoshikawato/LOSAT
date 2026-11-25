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
fn database_lookup_files_collects_alias_index_group_and_existing_dblist() {
    use super::database_lookup_files;
    // Plain nucleotide and protein names yield group alias + alias + index in lookup order.
    assert_eq!(
        database_lookup_files(b"x", SeqType::Nucleotide, None),
        vec![b"index.alx".to_vec(), b"x.nal".to_vec(), b"x.nin".to_vec()]
    );
    assert_eq!(
        database_lookup_files(b"p5", SeqType::Protein, None),
        vec![
            b"index.alx".to_vec(),
            b"p5.pal".to_vec(),
            b"p5.pin".to_vec()
        ]
    );
    // Multiple quoted names with spaces.
    assert_eq!(
        database_lookup_files(b"\"sub dir/x1\" x2", SeqType::Nucleotide, None),
        vec![
            b"sub dir/index.alx".to_vec(),
            b"sub dir/x1.nal".to_vec(),
            b"sub dir/x1.nin".to_vec(),
            b"index.alx".to_vec(),
            b"x2.nal".to_vec(),
            b"x2.nin".to_vec(),
        ]
    );
    // Existing alias file expands DBLIST and stops on mutual recursion.
    let dir = crate::utils::test_dir::TempDir::unique("losat-d19b-alias");
    std::fs::write(
        dir.join("a.nal"),
        b"TITLE a\nDBLIST b v1\nSEQIDLIST ids.txt\n",
    )
    .unwrap();
    std::fs::write(dir.join("b.nal"), b"TITLE b\nDBLIST a\n").unwrap();
    let base = dir.join("a");
    let got = database_lookup_files(
        base.as_os_str().as_encoded_bytes(),
        SeqType::Nucleotide,
        None,
    );
    let has = |name: &str| {
        let full = dir.join(name);
        got.iter()
            .any(|p| p.as_slice() == full.as_os_str().as_encoded_bytes())
    };
    assert!(has("a.nal"));
    assert!(has("b.nal"));
    assert!(has("v1.nal"));
    assert!(has("v1.nin"));
    assert!(has("ids.txt"));
    assert!(has("index.alx"));
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
