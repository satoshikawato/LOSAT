use super::*;

fn filter_of(num_oids: u32, included: &[u32]) -> OidFilter {
    let mut words = vec![0u64; (num_oids as usize).div_ceil(64)];
    for &oid in included {
        words[(oid / 64) as usize] |= 1 << (oid % 64);
    }
    OidFilter::new(words, num_oids)
}

/// A deterministic pseudo-random sequence (xorshift).
fn next(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

#[test]
fn rank_select_and_runs_agree_with_a_list() {
    let mut state = 0x2026_1030_u64;
    for &(num_oids, per_mille) in &[
        (1u32, 1000u64),
        (63, 500),
        (64, 1000),
        (65, 0),
        (511, 900),
        (512, 10),
        (513, 500),
        (5000, 2),
        (5000, 998),
        (20_000, 300),
    ] {
        let included: Vec<u32> = (0..num_oids)
            .filter(|_| next(&mut state) % 1000 < per_mille)
            .collect();
        let filter = filter_of(num_oids, &included);
        assert_eq!(filter.count() as usize, included.len(), "{num_oids}");
        for (ordinal, &oid) in included.iter().enumerate() {
            assert_eq!(
                filter.select(ordinal as u32),
                oid,
                "{num_oids} select {ordinal}"
            );
        }
        for oid in 0..=num_oids {
            let below = included.partition_point(|&o| o < oid) as u32;
            assert_eq!(filter.rank(oid), below, "{num_oids} rank {oid}");
            assert_eq!(filter.contains(oid), included.binary_search(&oid).is_ok());
            let limit = num_oids.min(oid + 130);
            let mut end = oid;
            while end < limit && included.binary_search(&end).is_ok() {
                end += 1;
            }
            assert_eq!(
                filter.run_end(oid, limit),
                end.max(oid.min(limit)),
                "{num_oids} run {oid}"
            );
        }
    }
}

/// NCBI's `x_CheckOrFindOID` (c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:272-296, cited at `OidFilter::check_or_find`)
/// over the list `listed` of `num_oids` OIDs with the range `[begin, end)` once set up.
fn ncbi_check_or_find(
    listed: &[u32],
    num_oids: u32,
    (begin, end): (u32, u32),
    set_up: &mut bool,
    next: &mut u32,
) -> bool {
    let (restrict_begin, restrict_end) = if *set_up { (begin, end) } else { (0, num_oids) };
    let mut success = true;
    if *next < restrict_begin {
        *next = restrict_begin;
    }
    if *next >= restrict_end {
        success = false;
    }
    *set_up = true;
    if success {
        match listed.iter().find(|&&oid| oid >= *next) {
            Some(&oid) => {
                *next = oid;
                success = oid <= end;
            }
            None => success = false,
        }
    }
    success
}

#[test]
fn the_range_and_ncbis_checks_agree_with_ncbis_loop() {
    let mut state = 0x2026_1113_u64;
    for round in 0..400u32 {
        let num_oids = [0u32, 1, 5, 64, 65, 130, 700][(round % 7) as usize];
        let per_mille = [0u64, 300, 700, 1000][(round / 7 % 4) as usize];
        let listed: Vec<u32> = (0..num_oids)
            .filter(|_| next(&mut state) % 1000 < per_mille)
            .collect();
        let end = (next(&mut state) % (u64::from(num_oids) + 1)) as u32;
        let begin = (next(&mut state) % (u64::from(end) + 1)) as u32;
        let filter = filter_of(num_oids, &listed).with_range((begin, end));
        let included: Vec<u32> = listed
            .iter()
            .copied()
            .filter(|oid| (begin..end).contains(oid))
            .collect();
        let case = format!("{num_oids} {listed:?} [{begin}, {end})");
        // The search's view: the listed OIDs inside the range.
        assert_eq!(filter.count() as usize, included.len(), "{case}");
        for (ordinal, &oid) in included.iter().enumerate() {
            assert_eq!(filter.select(ordinal as u32), oid, "{case}");
        }
        for oid in 0..=num_oids {
            let below = included.partition_point(|&o| o < oid) as u32;
            assert_eq!(filter.rank(oid), below, "{case} rank {oid}");
            assert_eq!(filter.contains(oid), included.contains(&oid), "{case}");
            let limit = num_oids.min(oid + 70);
            let mut stop = oid;
            while stop < limit && included.contains(&stop) {
                stop += 1;
            }
            assert_eq!(
                filter.run_end(oid, limit),
                stop.max(oid.min(limit)),
                "{case}"
            );
        }
        // NCBI's loop from OID 0, with the list set up before it or not.
        for set_up in [false, true] {
            let mut flag = set_up;
            let mut want = Vec::new();
            let mut oid = 0;
            while ncbi_check_or_find(&listed, num_oids, (begin, end), &mut flag, &mut oid) {
                want.push(oid);
                oid += 1;
            }
            let walk: Vec<u32> = filter.walk(set_up).oids(&filter).collect();
            assert_eq!(walk, want, "{case} set up {set_up}");
            // Each check alone.
            for oid in 0..=num_oids {
                let (mut flag, mut at, mut got) = (set_up, oid, oid);
                let want = ncbi_check_or_find(&listed, num_oids, (begin, end), &mut flag, &mut at);
                assert_eq!(
                    (filter.check_or_find(&mut got, set_up), got),
                    (want, at),
                    "{case} set up {set_up} check {oid}"
                );
            }
        }
    }
}

#[test]
fn bits_past_the_last_oid_are_not_included() {
    let filter = OidFilter::new(vec![u64::MAX; 2], 70);
    assert_eq!(filter.count(), 70);
    assert!(!filter.contains(70));
    assert_eq!(filter.run_end(60, 128), 70);
    assert_eq!(filter.select(69), 69);
}

fn tx_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/lmdb")
}

/// NCBI's OIDs of `tx` for each taxid set (`blastdbcmd -taxids S -no_taxid_expansion` for a
/// positive set, the sequences `blastp -negative_taxids S -no_taxid_expansion` searched for a
/// negative one: `tests/seqdb/lmdb/expected/tx_oracle.txt`).
#[test]
fn restricted_oids_are_ncbis_on_the_tx_fixture() {
    let oracle = std::fs::read_to_string(tx_dir().join("expected/tx_oracle.txt")).unwrap();
    let name = tx_dir().join("tx");
    let mut checked = 0;
    for line in oracle.lines().filter(|line| !line.trim().is_empty()) {
        let mut fields = line.split(' ').filter(|field| !field.is_empty());
        let kind = fields.next().unwrap();
        let set: BTreeSet<i32> = fields
            .next()
            .unwrap()
            .split(',')
            .map(|taxid| taxid.parse().unwrap())
            .collect();
        let rest: Vec<&str> = fields.collect();
        let list = match kind {
            "P" => TaxidList::Positive(set),
            "N" => TaxidList::Negative(set),
            other => panic!("{other}"),
        };
        let mut db = SeqDb::open(name.to_str().unwrap(), SeqType::Protein).unwrap();
        let whole = (db.num_seqs() as i32, db.total_length());
        let result = db.restrict_taxids(&list);
        if rest.first() == Some(&"Error:") {
            let Err(SeqDbError::NotFound(text)) = result else {
                panic!("{line}: {result:?}");
            };
            assert_eq!(
                text, b"Taxonomy ID(s) not found in the tx database.",
                "{line}"
            );
            continue;
        }
        let warning = result.unwrap_or_else(|error| panic!("{line}: {error}"));
        let want: Vec<u32> = rest.iter().map(|oid| oid.parse().unwrap()).collect();
        let got: Vec<u32> = (0..db.num_seqs())
            .map(|ordinal| db.oid(ordinal).get())
            .collect();
        assert_eq!(got, want, "{line}");
        assert_eq!(
            (db.printed_num_seqs(), db.printed_total_length()),
            whole,
            "{line}"
        );
        // A negative list of taxids that `tx` does not hold excludes nothing and warns.
        assert_eq!(
            warning.is_some(),
            matches!(&list, TaxidList::Negative(set) if set.iter().all(|t| *t == 999)),
            "{line}"
        );
        // The sequences read through the ordinal ids are those of the OIDs.
        let mut whole_db = SeqDb::open(name.to_str().unwrap(), SeqType::Protein).unwrap();
        whole_db.filter = None;
        let mut lengths = 0u64;
        for (ordinal, &oid) in want.iter().enumerate() {
            let filtered = db.read_sequence(ordinal as u32).unwrap();
            let plain = whole_db.read_sequence(oid).unwrap();
            assert_eq!(filtered.protein(0), plain.protein(0), "{line} {oid}");
            lengths += filtered.protein(0).len() as u64;
        }
        assert_eq!(db.total_length(), lengths, "{line}");
        let mut first = 0;
        let mut seen = Vec::new();
        let mut listed = Vec::new();
        while let Some(range) = db.next_block(first, 1 << 20) {
            listed.push(range);
            for at in range.first..range.end {
                seen.push(db.volume_starts[range.volume] + at);
            }
            first = range.end_ordinal();
        }
        assert_eq!(seen, want, "{line}: blocks");
        // LOSAT_X_BNDBSMALL: listed in runs, the blocks are the same (cut at excluded OIDs).
        assert_eq!(db.x_next_blocks(0, 1 << 20, usize::MAX), listed, "{line}");
        for max_bytes in [0, 1, 100] {
            let mut first = 0;
            let mut listed = Vec::new();
            while let Some(range) = db.next_block(first, max_bytes) {
                listed.push(range);
                first = range.end_ordinal();
            }
            for entries in [1, 2, 1 << 14] {
                assert_eq!(
                    db.x_next_blocks_in(0, max_bytes, usize::MAX, entries),
                    listed,
                    "{line} within {max_bytes}, runs of {entries}"
                );
            }
        }
        checked += 1;
    }
    assert!(checked >= 30, "{checked}");
}

fn ids_db(name: &str, seqtype: SeqType) -> SeqDb {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/db_search/ids")
        .join(name);
    SeqDb::open(path.to_str().unwrap(), seqtype).unwrap()
}

fn list_of(text: &[u8], negative: bool) -> SeqidList {
    let mut warnings = Vec::new();
    SeqidList {
        negative,
        file: seqidlist::read_si_list(text, &mut warnings).unwrap(),
    }
}

/// The Seq-id strings (`GetSeqIdString` with true) of the definition lines the restricted
/// database keeps, per ordinal id.
fn kept_ids(db: &SeqDb) -> Vec<Vec<String>> {
    (0..db.num_seqs())
        .map(|ordinal| {
            let at = db.locate_oid(db.oid(ordinal).get()).unwrap();
            let set = db.volumes[at.volume].defline_set(at.vol_oid).unwrap();
            set.deflines
                .iter()
                .flat_map(|d| d.ids.iter())
                .map(|id| {
                    let mut out = Vec::new();
                    label::content(id, true, &mut out);
                    String::from_utf8(out).unwrap()
                })
                .collect()
        })
        .collect()
}

#[test]
fn the_seq_align_id_lookup_finds_the_sequence_ncbi_finds_first() {
    // tests/db_search/ids/make_ids.py (line DBID N1, probe R/accoid/p1 with NCBI 2.17.0):
    // the Seq-align id of the second record of pairs 1 to 3 finds the first record, those of
    // pairs 4 and 5 find their own records; version 5 through the LMDB key, version 4 the
    // ISAM.
    // `pdup_ab` and `pdup_ac` are aliases of two databases (two LMDB files): the second's
    // record finds the first's.
    let pairs = [0, 0, 2, 2, 4, 4, 6, 7, 8, 9].map(Some);
    for (name, seqtype, want) in [
        ("pfound_v5", SeqType::Protein, &pairs[..]),
        ("pfound_v4", SeqType::Protein, &pairs[..]),
        ("nfound_v5", SeqType::Nucleotide, &pairs[..2]),
        ("nfound_v4", SeqType::Nucleotide, &pairs[..2]),
        ("pdup_ab", SeqType::Protein, &pairs[..2]),
        ("pdup_ac", SeqType::Protein, &pairs[..2]),
    ] {
        let db = ids_db(name, seqtype);
        let lookup = db.seqid_lookup().unwrap();
        let found: Vec<Option<u32>> = (0..db.num_seqs())
            .map(|ordinal| {
                let set = db.defline_set(ordinal).unwrap();
                lookup.seqid_to_oid(set.alignment_id().unwrap()).unwrap()
            })
            .collect();
        assert_eq!(found, want, "{name}");
    }
    // A list that leaves the first record out: `gb|AAA00001|` finds its own record, as
    // NCBI's lookup through `x_CheckOrFindOID` (row db.id.blastp.found.neg.6; version 4
    // keeps every OID here and refuses such a list when the search reads the volume).
    let mut db = ids_db("pfound_v5", SeqType::Protein);
    db.restrict(
        &DbList::Seqids(list_of(b"AAA00001.1\n", true)),
        &mut Vec::new(),
    )
    .unwrap();
    assert_eq!((db.oid(0).get(), db.ordinal_of(1)), (1, 0));
    let lookup = db.seqid_lookup().unwrap();
    let set = db.defline_set(0).unwrap();
    assert_eq!(
        lookup.seqid_to_oid(set.alignment_id().unwrap()).unwrap(),
        Some(1)
    );
}

#[test]
fn a_positive_seqid_list_keeps_the_named_oids_and_their_lines() {
    let whole = ids_db("prot_v5", SeqType::Protein);
    let mut db = ids_db("prot_v5", SeqType::Protein);
    let mut warnings = Vec::new();
    db.restrict(
        &DbList::Seqids(list_of(b">DAA10001.1\n", false)),
        &mut warnings,
    )
    .unwrap();
    assert_eq!(
        warnings,
        [ListWarning::Posted(TEXT_LIST_WARNING.to_string())]
    );
    assert!(db.is_restricted());
    assert_eq!(db.num_seqs(), 1);
    assert_eq!(kept_ids(&db), [vec!["DAA10001.1".to_string()]]);
    assert!(db.total_length() > 0 && db.total_length() < whole.total_length());
    assert_eq!(db.printed_num_seqs(), whole.num_seqs() as i32);

    // Versionless accessions, a pdb chain, a general id as db:tag, an unknown id, a gi.
    let mut db = ids_db("prot_v5", SeqType::Protein);
    db.restrict(
        &DbList::Seqids(list_of(
            b"NP_100006 XP_000011\n2XYZ_B mydb:prot10\nUNKNOWN1\ngi|1414\n",
            false,
        )),
        &mut Vec::new(),
    )
    .unwrap();
    assert_eq!(db.num_seqs(), 4);
    for ids in kept_ids(&db) {
        assert!(
            ids.iter().any(|id| {
                ["NP_100006.1", "XP_000011.2", "2XYZ_B", "mydb:prot10"].contains(&id.as_str())
            }),
            "{ids:?}"
        );
    }
}

#[test]
fn a_negative_seqid_list_excludes_an_oid_only_when_every_id_is_covered() {
    let whole = ids_db("prot_v5", SeqType::Protein).num_seqs();
    // `lcl|Subject_28|tpg|DAA10008.1|`: one id listed keeps the OID, both exclude it (the
    // versionless `DAA10008` covers the stored `DAA10008.1`).
    let mut db = ids_db("prot_v5", SeqType::Protein);
    db.restrict(
        &DbList::Seqids(list_of(b"Subject_28\n", true)),
        &mut Vec::new(),
    )
    .unwrap();
    assert_eq!(db.num_seqs(), whole);
    assert!(db.is_restricted());
    let mut db = ids_db("prot_v5", SeqType::Protein);
    db.restrict(
        &DbList::Seqids(list_of(b"Subject_28\nDAA10008\n", true)),
        &mut Vec::new(),
    )
    .unwrap();
    assert_eq!(db.num_seqs(), whole - 1);
    // A definition line whose only ids are gis or `ti` general ids is dropped by a
    // negative seqid list (`CSeqDBNegativeList::FindId` matches no type): an OID of such
    // lines alone keeps none: its hits get no Seq-align (`keeps_line`), and a read of its
    // lines still refuses (no report reads it); every other kept OID reads.
    let unrestricted = ids_db("prot_v5", SeqType::Protein);
    let gi_only: Vec<u32> = (0..unrestricted.num_seqs())
        .filter(|&oid| {
            let at = unrestricted.locate_oid(oid).unwrap();
            let set = unrestricted.volumes[at.volume]
                .defline_set(at.vol_oid)
                .unwrap();
            set.deflines.iter().all(|d| {
                d.ids.iter().all(|id| {
                    matches!(id, ParsedSeqId::Gi(_))
                        || matches!(id, ParsedSeqId::General(tag) if tag.db == b"ti")
                })
            })
        })
        .collect();
    assert!(!gi_only.is_empty());
    for ordinal in 0..db.num_seqs() {
        let oid = db.oid(ordinal).get();
        let at = db.locate_oid(oid).unwrap();
        match db.volumes[at.volume].defline_set(at.vol_oid) {
            Ok(set) => {
                assert!(!gi_only.contains(&oid), "{oid}");
                assert!(db.keeps_line(ordinal).unwrap(), "{oid}");
                assert!(set.deflines.iter().all(|d| {
                    d.ids.iter().all(|id| {
                        let mut out = Vec::new();
                        label::content(id, true, &mut out);
                        out != b"DAA10008.1"
                    })
                }));
            }
            Err(error) => {
                assert!(gi_only.contains(&oid), "{oid}: {error}");
                assert!(error.to_string().contains("none of its definition lines"));
                assert!(!db.keeps_line(ordinal).unwrap(), "{oid}");
            }
        }
    }
    // A nucleotide local id of its own OID.
    let whole = ids_db("nucl_v5", SeqType::Nucleotide).num_seqs();
    let mut db = ids_db("nucl_v5", SeqType::Nucleotide);
    db.restrict(
        &DbList::Seqids(list_of(b"mynuc6\nUNKNOWN1\n", true)),
        &mut Vec::new(),
    )
    .unwrap();
    assert_eq!(db.num_seqs(), whole - 1);
}

#[test]
fn a_list_that_loses_every_id_leaves_a_version_5_database_of_one_volume_whole() {
    let whole = ids_db("prot_v5", SeqType::Protein).num_seqs();
    for negative in [false, true] {
        let mut db = ids_db("prot_v5", SeqType::Protein);
        let mut warnings = Vec::new();
        db.restrict(
            &DbList::Seqids(list_of(b"gi|1414\ngi|2\n", negative)),
            &mut warnings,
        )
        .unwrap();
        assert_eq!(
            warnings,
            [
                ListWarning::Posted(TEXT_LIST_WARNING.to_string()),
                ListWarning::Posted(EMPTY_LIST_WARNING.to_string())
            ]
        );
        assert!(db.is_restricted());
        assert_eq!(db.num_seqs(), whole);
        assert_eq!(kept_ids(&db).len(), whole as usize);
        assert!(db.min_length() > 0);
    }
}

// Line DBID gisl: `x_ApplyUserGiList` clears every OID of a positive list that lost every id,
// except where NCBI's OID bit set is the special all-set (one volume that the alias files
// keep whole); a negative one excludes nothing.
#[test]
fn a_positive_list_that_loses_every_id_keeps_no_oid_unless_one_volume_is_kept_whole() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    for (name, seqtype, whole) in [
        ("taxfilter/tfp", SeqType::Protein, true),
        // `FIRST_OID 1`, `LAST_OID 23`: the whole volume.
        ("aliaskeys/p_r123", SeqType::Protein, true),
        // A range and the volume itself; the same volume twice.
        ("aliaskeys/p_twice", SeqType::Protein, true),
        ("aliaskeys/p_dupvol", SeqType::Protein, true),
        // The whole volume and a range inside it, in both orders.
        ("aliaskeys/p_cover", SeqType::Protein, true),
        ("aliaskeys/p_coverr", SeqType::Protein, true),
        // An alias of an alias; the whole range above it.
        ("aliaskeys/p_chain", SeqType::Protein, true),
        ("aliaskeys/p_above", SeqType::Protein, true),
        // Two ranges ORed (a bit array), then the whole range.
        ("aliaskeys/p_arrcov", SeqType::Protein, true),
        // One volume of a nucleotide database.
        ("taxfilter/tfn.00", SeqType::Nucleotide, true),
        // A range inside the volume; two ranges ORed, and two whose union is the volume.
        ("aliaskeys/p_range", SeqType::Protein, false),
        ("aliaskeys/p_union", SeqType::Protein, false),
        ("aliaskeys/p_halves", SeqType::Protein, false),
        // A range that cuts the volume above an alias without one.
        ("aliaskeys/p_abovecut", SeqType::Protein, false),
        // Two volumes; one of them inside a range.
        ("taxfilter/tfn", SeqType::Nucleotide, false),
        ("aliaskeys/n_vol1", SeqType::Nucleotide, false),
    ] {
        let path = root.join(name);
        let path = path.to_str().unwrap();
        let unrestricted = SeqDb::open(path, seqtype).unwrap().num_seqs();
        for negative in [false, true] {
            let mut db = SeqDb::open(path, seqtype).unwrap();
            let mut warnings = Vec::new();
            db.restrict(
                &DbList::Seqids(list_of(b"gi|1414\n", negative)),
                &mut warnings,
            )
            .unwrap();
            assert_eq!(
                warnings.last(),
                Some(&ListWarning::Posted(EMPTY_LIST_WARNING.to_string())),
                "{name}"
            );
            let kept = if negative || whole { unrestricted } else { 0 };
            assert_eq!(db.num_seqs(), kept, "{name}, negative {negative}");
        }
    }
}

// Probes T08-T12, C01-C97 of line DI-8: 0 OIDs, NCBI's statistics read OID 0 and one sequence.
#[test]
fn a_list_without_a_kept_sequence_has_the_statistics_of_oid_zero() {
    for (name, seqtype) in [
        ("prot_v5", SeqType::Protein),
        ("nucl_v5", SeqType::Nucleotide),
    ] {
        let mut db = ids_db(name, seqtype);
        let first = match seqtype {
            SeqType::Protein => u64::from(db.protein_length(0).unwrap()),
            SeqType::Nucleotide => db.nucleotide_length(0).unwrap() as u64,
        };
        let whole = (db.total_length() as i64, db.num_seqs() as i32);
        assert!(!db.filtered_empty());
        assert_eq!(db.statistics_numbers(), whole, "{name}");
        let mut warnings = Vec::new();
        db.restrict(
            &DbList::Seqids(list_of(b"UNKNOWN1\n", false)),
            &mut warnings,
        )
        .unwrap();
        assert!(db.filtered_empty(), "{name}");
        assert_eq!((db.num_seqs(), db.total_length()), (0, 0), "{name}");
        assert_eq!(db.statistics_numbers(), (first as i64, 1), "{name}");
        assert_eq!(db.printed_num_seqs(), whole.1, "{name}");
        assert!(db.next_block(0, 1 << 20).is_none(), "{name}");
    }
}

#[test]
fn a_binary_list_selects_what_its_text_list_selects() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/seqidlist");
    let text = std::fs::read(dir.join("sl_nucl.txt")).unwrap();
    let binary = std::fs::read(dir.join("sl_nucl.bsl")).unwrap();
    let mut from_text = ids_db("nucl_v5", SeqType::Nucleotide);
    let mut warnings = Vec::new();
    from_text
        .restrict(&DbList::Seqids(list_of(&text, false)), &mut warnings)
        .unwrap();
    assert_eq!(warnings.len(), 1);
    let mut from_binary = ids_db("nucl_v5", SeqType::Nucleotide);
    let mut warnings = Vec::new();
    from_binary
        .restrict(&DbList::Seqids(list_of(&binary, false)), &mut warnings)
        .unwrap();
    assert!(warnings.is_empty());
    assert!(from_text.num_seqs() >= 6);
    assert_eq!(from_text.num_seqs(), from_binary.num_seqs());
    assert_eq!(kept_ids(&from_text), kept_ids(&from_binary));
    assert_eq!(from_text.total_length(), from_binary.total_length());
    // A recorded database length that is not this database's: NCBI's warning, the list used.
    let badlen = std::fs::read(dir.join("sl_tfp_badlen.bsl")).unwrap();
    let mut db = ids_db("prot_v5", SeqType::Protein);
    let mut warnings = Vec::new();
    let result = db.restrict(&DbList::Seqids(list_of(&badlen, false)), &mut warnings);
    assert_eq!(warnings, [ListWarning::Posted(DB_INFO_WARNING.to_string())]);
    // tfp's ids are not in prot_v5: no OID is kept (the search programs refuse that).
    assert!(result.is_ok() && db.num_seqs() == 0);
}

#[test]
fn a_sequence_kept_without_a_definition_line_gets_no_seq_align() {
    // `lcl|Subject_29|pdb|3ABC|C`: the molecule name is a lookup key, no line's id. The
    // list keeps the sequence (it is searched and counted), but none of its lines: NCBI
    // makes no Seq-align of its hits (`keeps_line`); a read of its lines refuses.
    let mut db = ids_db("prot_v5", SeqType::Protein);
    db.restrict(
        &DbList::Seqids(list_of(b"3ABC\nDAA10001.1\n", false)),
        &mut Vec::new(),
    )
    .unwrap();
    assert_eq!(db.num_seqs(), 2);
    assert!(db.drops_lineless_hits());
    let kept: Vec<bool> = (0..2)
        .map(|ordinal| db.keeps_line(ordinal).unwrap())
        .collect();
    assert_eq!(kept.iter().filter(|&&kept| kept).count(), 1, "{kept:?}");
    let lineless = kept.iter().position(|&kept| !kept).unwrap() as u32;
    let error = db.defline_set(lineless).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("but none of its definition lines, and NCBI leaves its hits out"),
        "{error}"
    );
    assert!(error.ncbi_message().is_none());
    // Without a seqid list every sequence keeps its lines, and none is read.
    let whole = ids_db("prot_v5", SeqType::Protein);
    assert!(!whole.drops_lineless_hits());
    assert!(whole.keeps_line(0).unwrap());
}

#[test]
fn a_negative_list_with_ids_excludes_its_lineless_sequences_at_the_fetches() {
    // A negative seqid list keeps a line with an id that is neither a gi nor a ti nor listed.
    // prot_v4's sequences 14 (`gi|1414`) and 28 (`gnl|ti|123`) have one line, with a gi or a
    // ti only, and the one line of sequence 11, `ref|XP_000011.2|`, is the list's
    // `XP_000011`: the list keeps the three sequences (NCBI 2.17.0's statistics count all 53
    // under it) and removes their lines. NCBI's id set holds the list's ids, so its fetches
    // with `check_oid_exclusion` exclude them.
    let mut db = ids_db("prot_v4", SeqType::Protein);
    db.restrict(
        &DbList::Seqids(list_of(b"XP_000011\n", true)),
        &mut Vec::new(),
    )
    .unwrap();
    assert!(db.drops_lineless_hits());
    assert!(db.excludes_lineless_at_fetch());
    let lineless: Vec<u32> = (0..db.num_seqs() as u32)
        .filter(|&ordinal| !db.keeps_line(ordinal).unwrap())
        .collect();
    assert_eq!(lineless, [11, 14, 28]);
    // A positive list leaves NCBI's id set blank: its sequences without a line are not
    // excluded at the fetches, and leave at the Seq-aligns.
    let mut positive = ids_db("prot_v5", SeqType::Protein);
    positive
        .restrict(
            &DbList::Seqids(list_of(b"3ABC\nDAA10001.1\n", false)),
            &mut Vec::new(),
        )
        .unwrap();
    assert!(positive.drops_lineless_hits());
    assert!(!positive.excludes_lineless_at_fetch());
    // A negative list of gis only is empty on a version 5 database once the gis are dropped
    // (`s_ProcessSeqIDsForV5`): no line rule, no id set.
    let mut gis = ids_db("prot_v5", SeqType::Protein);
    gis.restrict(
        &DbList::Seqids(list_of(b"gi|1414\n", true)),
        &mut Vec::new(),
    )
    .unwrap();
    assert!(!gis.drops_lineless_hits());
    assert!(!gis.excludes_lineless_at_fetch());
}

#[test]
fn a_version_4_database_refuses_a_binary_seqid_list() {
    let binary = std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/seqidlist/sl_tfp.bsl"),
    )
    .unwrap();
    let mut db = ids_db("prot_v4", SeqType::Protein);
    let error = db
        .restrict(&DbList::Seqids(list_of(&binary, false)), &mut Vec::new())
        .unwrap_err();
    assert_eq!(
        error.ncbi_message().as_deref(),
        Some(BINARY_LIST_V4_ERROR.as_bytes())
    );
}

#[test]
fn a_version_4_database_refuses_a_taxid_list() {
    let name = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db/prot_v4");
    let mut db = SeqDb::open(name.to_str().unwrap(), SeqType::Protein).unwrap();
    let error = db
        .restrict_taxids(&TaxidList::Positive([9606].into_iter().collect()))
        .unwrap_err();
    assert_eq!(
        error.ncbi_message().as_deref(),
        Some(b"Taxonomy filtering is not supported in v4 BLAST dbs".as_slice())
    );
}

fn defline(taxid: Option<i32>, links: Option<Vec<i32>>) -> BlastDefline {
    BlastDefline {
        taxid,
        links,
        ..BlastDefline::default()
    }
}

#[test]
fn definition_lines_follow_ncbis_rules() {
    let set = |ids: &[i32]| ids.iter().copied().collect::<BTreeSet<i32>>();
    let positive = DeflineRule::Taxids(set(&[9606, 10090]));
    assert!(positive.keeps(&defline(Some(9606), None)));
    assert!(positive.keeps(&defline(Some(1), Some(vec![10090]))));
    assert!(!positive.keeps(&defline(Some(1), Some(vec![2]))));
    // `IsSetTaxid`: an absent taxid is not 0.
    assert!(!DeflineRule::Taxids(set(&[0])).keeps(&defline(None, None)));
    assert!(DeflineRule::Taxids(set(&[0])).keeps(&defline(Some(0), None)));

    let negative = DeflineRule::NegativeTaxids(set(&[9606]));
    assert!(!negative.keeps(&defline(Some(9606), None)));
    assert!(negative.keeps(&defline(Some(9606), Some(vec![10090]))));
    // `GetTaxIds` drops 0 beside another taxid.
    assert!(!negative.keeps(&defline(Some(0), Some(vec![9606]))));
    assert!(negative.keeps(&defline(Some(0), None)));
    // No taxid at all: an empty set is a subset of every list.
    assert!(!negative.keeps(&defline(None, None)));
}

#[test]
fn a_rule_that_keeps_no_line_keeps_every_line() {
    let rule = DeflineRule::Taxids([5].into_iter().collect());
    let set = DeflineSet {
        deflines: vec![defline(Some(1), None), defline(Some(5), None)],
        utf8: false,
    };
    assert_eq!(rule.keep_mask(&set), Some(vec![false, true]));
    let none = DeflineSet {
        deflines: vec![defline(Some(1), None), defline(Some(2), None)],
        utf8: false,
    };
    assert_eq!(rule.keep_mask(&none), Some(vec![true, true]));
    let mut lines = vec![1, 2, 3];
    lines.retain_lines(&[true, false, true]);
    assert_eq!(lines, vec![1, 3]);
    lines.retain_lines(&[false]);
    assert_eq!(lines, vec![1, 3]);
}

/// NCBI BLAST+ 2.17.0's texts of an LMDB file and of a lookup file that does not open or whose
/// bytes cannot be had (line WQA D18 probes 5-7, 2026-11-21): the LMDB error for every code,
/// EMFILE included, and `CMemoryFile`'s exception for the lookup files.
#[cfg(target_os = "linux")]
#[test]
fn unreadable_lmdb_and_lookup_files_give_ncbis_errors() {
    let open = |code: i32| SeqDbError::Open {
        path: PathBuf::from("/d/tm_psi.pos"),
        source: std::io::Error::from_raw_os_error(code),
    };
    // The first byte, read where the file is opened (`DbFile::open`): NCBI's mapping.
    let read = |code: i32| SeqDbError::Map {
        path: PathBuf::from("/d/tm_psi.pos"),
        len: 1234,
        source: std::io::Error::from_raw_os_error(code),
    };
    let lmdb = |error: SeqDbError| match lmdb_file_error(error) {
        SeqDbError::NotFound(text) => String::from_utf8(text).unwrap(),
        other => panic!("{other:?}"),
    };
    assert_eq!(
        lmdb(open(13)),
        "LMDB runtime error: mdb_env_open: Permission denied"
    );
    assert_eq!(
        lmdb(open(24)),
        "LMDB runtime error: mdb_env_open: Too many open files"
    );
    assert_eq!(
        lmdb(read(5)),
        "LMDB runtime error: mdb_env_open: Input/output error"
    );
    let lookup = |error: SeqDbError| match lookup_file_error(error) {
        SeqDbError::Ncbi(exception) => String::from_utf8(exception.what()).unwrap(),
        other => panic!("{other:?}"),
    };
    // NCBI's report in LOSAT's words (`crate::wording::losat_text`).
    let report = |line: u32, function: &str, text: &str| {
        crate::wording::losat_text(&format!(
            "NCBI C++ Exception:\n    T0 \"{}/src/corelib/ncbifile.cpp\", line {line}: Error: \
             (CFileException::eMemoryMap) ncbi::{function}() - {text}\n",
            crate::ncbi_exception::NCBI_BUILD_ROOT
        ))
        .into_owned()
    };
    assert_eq!(
        lookup(open(13)),
        report(
            5945,
            "CMemoryFileMap::x_Open",
            "CMemoryFile: Cannot memory map file '/d/tm_psi.pos':Permission denied"
        )
    );
    assert_eq!(
        lookup(read(5)),
        report(
            5650,
            "CMemoryFileSegment::CMemoryFileSegment",
            "Cannot map file '/d/tm_psi.pos' to memory (offset=0, length=1234): Input/output error"
        )
    );
    // A first byte that cannot be read without a system code, and a read that fails after the
    // first byte (where NCBI ends by SIGBUS), are LOSAT's own errors.
    let short = SeqDbError::Map {
        path: PathBuf::from("/d/tm_psi.pos"),
        len: 1234,
        source: std::io::Error::from(std::io::ErrorKind::UnexpectedEof),
    };
    assert!(matches!(lookup_file_error(short), SeqDbError::Map { .. }));
    let later = || SeqDbError::Read {
        path: PathBuf::from("/d/tm_psi.pos"),
        source: std::io::Error::from_raw_os_error(5),
    };
    assert!(matches!(
        lookup_file_error(later()),
        SeqDbError::Read { .. }
    ));
    assert!(matches!(lmdb_file_error(later()), SeqDbError::Read { .. }));
}

/// A read of an LMDB or lookup file that fails after the file opened (its first byte was read) is
/// LOSAT's own error, "<path> cannot be read" (`SeqDbError::Read`), not a refusal of the file's
/// format: NCBI maps these files, so such a failure is a SIGBUS there (line WQA D18c,
/// 2026-11-21). Each file of `tx` is read through a handle that holds only its first bytes
/// (`file::refused_reads`): the `.pdb` of a seqid list and the `.ptf` of a taxid list keep their
/// two meta pages (`Env::open` reads them; the open's volume list or the first lookup reads a page
/// past them), the `.pto` its
/// first byte, the `.pot` and `.pos` of the negative lists their OID count.
#[cfg(target_os = "linux")]
#[test]
fn a_later_read_of_an_lmdb_or_lookup_file_is_losats_own_error() {
    let page = Env::open(&std::fs::read(tx_dir().join("tx.pdb")).unwrap())
        .unwrap()
        .page_size();
    let taxids = |negative: bool| {
        let set: BTreeSet<i32> = [100].into_iter().collect();
        if negative {
            TaxidList::Negative(set)
        } else {
            TaxidList::Positive(set)
        }
    };
    let cases: [(&str, usize, Option<TaxidList>, &[u8], bool); 5] = [
        ("pdb", 2 * page, None, b"AAA00001.1\n", false),
        ("ptf", 2 * page, Some(taxids(false)), b"", false),
        ("pto", 1, Some(taxids(false)), b"", false),
        ("pot", 8, Some(taxids(true)), b"", false),
        ("pos", 8, None, b"AAA00001.1\n", true),
    ];
    for (ext, keep, taxid_list, seqids, negative) in cases {
        let dir = crate::utils::test_dir::temp_dir().join(format!(
            "losat-filter-later-read-{ext}-{}",
            crate::utils::test_dir::process_id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        for entry in std::fs::read_dir(tx_dir()).unwrap() {
            let path = entry.unwrap().path();
            let name = path.file_name().unwrap().to_str().unwrap().to_string();
            if name.starts_with("tx.") {
                std::fs::copy(&path, dir.join(&name)).unwrap();
            }
        }
        let file = dir.join(format!("tx.{ext}"));
        let short = dir.join("short");
        std::fs::write(&short, &std::fs::read(&file).unwrap()[..keep]).unwrap();
        crate::seqdb::file::refused_reads::refuse(&file, &short);
        // The open reads the volume list of the `.pdb` past its meta pages.
        let error = SeqDb::open(dir.join("tx").to_str().unwrap(), SeqType::Protein)
            .and_then(|mut db| match &taxid_list {
                Some(list) => db.restrict_taxids(list).map(drop),
                None => db
                    .restrict(&DbList::Seqids(list_of(seqids, negative)), &mut Vec::new())
                    .map(drop),
            })
            .expect_err(ext);
        assert!(
            matches!(&error, SeqDbError::Read { path, .. } if path == &file),
            "{ext}: {error:?}"
        );
        assert!(error.is_unreadable(), "{ext}");
        assert!(
            error
                .to_string()
                .starts_with(&format!("{} cannot be read: ", file.display())),
            "{ext}: {error}"
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}

/// NCBI's `mdb_env_open` reads the meta pages of an LMDB file with `pread` before it maps the
/// file: a read that fails there with a system code gives its text, as an open does; one without
/// a code is LOSAT's own error (line WQA D18c, 2026-11-21).
#[test]
fn a_meta_page_that_cannot_be_read_gives_mdb_env_opens_text() {
    let path = tx_dir().join("tx.pdb");
    let file = DbFile::open(&path).unwrap();
    let unreadable = |source: std::io::Error| {
        ReadError::Unreadable(SeqDbError::Read {
            path: path.clone(),
            source,
        })
    };
    match lmdb_open_error(&file, unreadable(std::io::Error::from_raw_os_error(5))) {
        SeqDbError::NotFound(text) => assert_eq!(
            text,
            b"LMDB runtime error: mdb_env_open: Input/output error"
        ),
        other => panic!("{other:?}"),
    }
    let error = lmdb_open_error(
        &file,
        unreadable(std::io::Error::from(std::io::ErrorKind::UnexpectedEof)),
    );
    assert!(error.is_unreadable(), "{error:?}");
    assert!(matches!(
        lmdb_open_error(&file, ReadError::Format("LMDB: x (MDB_INVALID)".into())),
        SeqDbError::NotFound(_)
    ));
}

#[test]
fn ordinals_in_names_the_ordinals_of_an_oid_range() {
    let included = [0u32, 3, 4, 9, 64, 65, 130];
    let filter = filter_of(140, &included);
    let map = OidMap(Some(&filter));
    for from in 0..=140 {
        for to in 0..=140 {
            let ordinals = map.ordinals_in(from, to);
            let expected: Vec<u32> = (0..included.len() as u32)
                .filter(|&ordinal| (from..to).contains(&included[ordinal as usize]))
                .collect();
            assert_eq!(ordinals.collect::<Vec<_>>(), expected, "{from}..{to}");
        }
    }
    assert_eq!(OidMap::unrestricted().ordinals_in(3, 7), 3..7);
    assert_eq!(OidMap::unrestricted().ordinals_in(7, 3), 7..7);
}
