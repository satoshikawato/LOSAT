use super::*;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
fn make_hsp(e_value: f64, raw_score: i32, s_idx: u32) -> BlastnHsp {
    BlastnHsp {
        identity: 0.0,
        length: 0,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 1,
        s_start: 1,
        s_end: 1,
        e_value,
        bit_score: 0.0,
        num_ident: 0,
        query_frame: 1,
        query_length: 1,
        q_idx: 0,
        s_idx,
        raw_score,
        internal_q_offset_0: 0,
        internal_q_end_0: 1,
        internal_s_offset_0: 0,
        internal_s_end_0: 1,
        internal_query_context_offset: 0,
        gap_info: None,
        num_positives: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1355
#[test]
fn stable_hsp_sort_preserves_comparator_equal_input_order() {
    let mut hsps = vec![
        make_hsp(1.0, 100, 7),
        make_hsp(1.0, 100, 3),
        make_hsp(1.0, 100, 9),
    ];

    qsort_blastn_hsps_by(&mut hsps, |_, _| Ordering::Equal);

    assert_eq!(
        hsps.iter().map(|hsp| hsp.s_idx).collect::<Vec<_>>(),
        vec![7, 3, 9]
    );
}

#[test]
fn test_prelim_hitlist_size_gapped() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:43-70
    assert_eq!(get_prelim_hitlist_size(1, false, true), 10);
    assert_eq!(get_prelim_hitlist_size(30, false, true), 60);
    assert_eq!(get_prelim_hitlist_size(1000, false, true), 1050);
    assert_eq!(
        get_prelim_hitlist_size((1 << 30) - 1, false, true),
        (1 << 30) + 49
    );
    // NCBI's `Int4` overflows from 2^30 on (NCBI BLAST+ 2.17.0: 10, then a negative size);
    // LOSAT keeps the mathematical size (PD-LOSAT-NCBI-DEFECTS, hit list sizes past `Int4`).
    assert_eq!(
        get_prelim_hitlist_size(1 << 30, false, true),
        (1 << 30) + 50
    );
    assert_eq!(
        get_prelim_hitlist_size(i32::MAX as usize, false, true) as u64,
        i32::MAX as u64 + 50
    );
}

#[test]
fn prelim_hitlist_size_with_composition_based_statistics() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:56-63
    assert_eq!(get_prelim_hitlist_size(1, true, true), 1050);
    assert_eq!(get_prelim_hitlist_size(500, true, true), 1050);
    assert_eq!(get_prelim_hitlist_size(501, true, true), 1052);
    assert_eq!(
        get_prelim_hitlist_size(1_073_741_798, true, true),
        2_147_483_646
    );
    // NCBI's `2 * h + 50` overflows its `Int4` from 1073741799 on (NCBI BLAST+ 2.17.0: a size of
    // 0 or less, then 2 to 48 from 2147483624); LOSAT keeps the mathematical size.
    assert_eq!(
        get_prelim_hitlist_size(1_073_741_799, true, true) as u64,
        2_147_483_648
    );
    assert_eq!(
        get_prelim_hitlist_size(i32::MAX as usize, true, true) as u64,
        4_294_967_344
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/blasthits_unit_test.cpp:1258-1292 testSBlastHitsParamsUngapped
// NCBI's case is a blastp search with `gapped_calculation` false and the default hit list
// size; the function is the same one (`GetPrelimHitlistSize`), which keeps the hit list size
// for an ungapped search.
#[test]
fn ungapped_prelim_hitlist_size_is_the_hit_list_size() {
    assert_eq!(get_prelim_hitlist_size(500, false, false), 500);
    assert_eq!(get_prelim_hitlist_size(1, false, false), 1);
    assert_eq!(
        get_prelim_hitlist_size(i32::MAX as usize, false, false),
        i32::MAX as usize
    );
    // A gapped search adds 50 (testSBlastHitsParamsGapped).
    assert_eq!(get_prelim_hitlist_size(500, false, true), 550);
}

#[test]
fn test_hitlist_update_keeps_best() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3297
    let mut hit_list = BlastnHitList::new(1);

    let hsp_list_a = BlastnHspList {
        oid: 10,
        query_index: 0,
        hsps: vec![make_hsp(5.0, 50, 10)],
        best_evalue: i32::MAX as f64,
    };
    hit_list.update(hsp_list_a);
    assert_eq!(hit_list.hsplist_count, 1);
    assert_eq!(hit_list.hsplist_array[0].oid, 10);

    let hsp_list_b = BlastnHspList {
        oid: 20,
        query_index: 0,
        hsps: vec![make_hsp(1.0, 80, 20)],
        best_evalue: i32::MAX as f64,
    };
    hit_list.update(hsp_list_b);
    assert_eq!(hit_list.hsplist_count, 1);
    assert_eq!(hit_list.hsplist_array[0].oid, 20);

    let hsp_list_c = BlastnHspList {
        oid: 30,
        query_index: 0,
        hsps: vec![make_hsp(10.0, 40, 30)],
        best_evalue: i32::MAX as f64,
    };
    hit_list.update(hsp_list_c);
    assert_eq!(hit_list.hsplist_count, 1);
    assert_eq!(hit_list.hsplist_array[0].oid, 20);
}

fn make_list(oid: u32, hsps: &[(f64, i32)]) -> BlastnHspList {
    BlastnHspList {
        oid,
        query_index: 0,
        hsps: hsps
            .iter()
            .map(|&(e_value, score)| make_hsp(e_value, score, oid))
            .collect(),
        best_evalue: i32::MAX as f64,
    }
}

fn list_oids(hit_list: &BlastnHitList) -> Vec<u32> {
    hit_list.hsplist_array.iter().map(|list| list.oid).collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3095-3106
// An equal list replaces the heap root (`evalue_order < 0` frees only a worse one,
// blast_hits.c:3285-3291), so equal lists keep the higher ordinal ids, first.
#[test]
fn hit_list_equal_evalue_and_score_keeps_higher_oids() {
    let mut hit_list = BlastnHitList::new(3);
    for oid in 0..10 {
        hit_list.update(make_list(oid, &[(1e-10, 100)]));
    }
    assert!(hit_list.heapified);
    hit_list.sort_by_evalue();
    assert_eq!(list_oids(&hit_list), vec![9, 8, 7]);
}

// NCBI 2.17.0 blastp -comp_based_stats 0 over `LENGTH 0` (probe p5_L0#57 of line RG D1a-3):
// without the final sort, a full hit list keeps the heap's order of the lists inserted by OID.
// Here: lists come in any order; the traceback order inserts them by OID into a hit list of
// 3 (the fourth list heapifies it, worst first), and drops the lists without HSPs.
#[test]
fn traceback_order_inserts_by_oid_then_keeps_the_heap() {
    let mut hit_list = BlastnHitList::new(10);
    for (oid, score) in [(4, 40), (1, 10), (3, 30), (0, 50), (2, 20)] {
        hit_list.update(make_list(oid, &[(10f64.powi(-score), score)]));
    }
    hit_list.update(BlastnHspList {
        oid: 5,
        query_index: 0,
        hsps: Vec::new(),
        best_evalue: 0.0,
    });
    let mut fits = hit_list.clone();
    fits.put_in_traceback_order(10);
    assert_eq!(list_oids(&fits), vec![0, 1, 2, 3, 4]);
    hit_list.put_in_traceback_order(3);
    assert!(hit_list.heapified);
    // OIDs 0, 1, 2 come first; 3 replaces the root (OID 1, score 10), the heap of 0, 3, 2
    // has 2 (score 20) at its root; 4 replaces it.
    let mut reference = BlastnHitList::new(3);
    for oid in 0..5 {
        let score = [50, 10, 20, 30, 40][oid as usize];
        reference.update(make_list(oid, &[(10f64.powi(-score), score)]));
    }
    assert_eq!(list_oids(&hit_list), list_oids(&reference));
    assert_eq!(hit_list.hsplist_array[0].oid, 3);
}

// LOSAT_X_BNHLGROW (on unless `0`): the array doubles from 4 as lists come and never holds
// more room than NCBI's array (`hsplist_current`: 100, 200, then `hsplist_max`); without the
// switch it is reserved to that size when a list comes.
#[test]
fn hit_list_array_grows_as_lists_come_within_ncbi_array_size() {
    let grows = x_bnhlgrow();
    let mut hit_list = BlastnHitList::new(250);
    let mut rooms = Vec::new();
    for oid in 0..600 {
        hit_list.update(make_list(oid, &[(1e-10, 100 + (oid % 13) as i32)]));
        let room = hit_list.hsplist_array.capacity();
        assert!(room >= hit_list.hsplist_count, "list {oid}");
        if grows {
            assert!(room <= hit_list.hsplist_current, "list {oid}: {room}");
        }
        rooms.push(room);
    }
    assert!(hit_list.heapified);
    assert_eq!(hit_list.hsplist_count, 250);
    if grows {
        assert_eq!(
            [rooms[0], rooms[3], rooms[4], rooms[64], rooms[99], rooms[100], rooms[200]],
            [4, 4, 8, 100, 100, 200, 250]
        );
    } else {
        assert_eq!([rooms[0], rooms[100]], [100, 200]);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1392-1395
#[test]
fn hit_list_evalues_below_1e180_compare_equal() {
    let mut hit_list = BlastnHitList::new(1);
    hit_list.update(make_list(1, &[(1e-200, 50)]));
    hit_list.update(make_list(2, &[(1e-190, 60)]));
    assert_eq!(list_oids(&hit_list), vec![2]);
    let mut hit_list = BlastnHitList::new(1);
    hit_list.update(make_list(1, &[(1e-100, 50)]));
    hit_list.update(make_list(2, &[(1e-90, 60)]));
    assert_eq!(list_oids(&hit_list), vec![1]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3099-3102
#[test]
fn hit_list_equal_evalue_higher_first_score_wins() {
    let mut hit_list = BlastnHitList::new(1);
    hit_list.update(make_list(5, &[(1e-5, 40)]));
    hit_list.update(make_list(1, &[(1e-5, 41)]));
    assert_eq!(list_oids(&hit_list), vec![1]);
    hit_list.update(make_list(9, &[(1e-5, 40)]));
    assert_eq!(list_oids(&hit_list), vec![1]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3276-3287
#[test]
fn hit_list_heap_keeps_the_best_lists_in_comparator_order() {
    let keys: Vec<(u32, f64, i32)> = (0..40u32)
        .map(|oid| {
            let e_value = [1e-30, 1e-10, 1e-3, 0.5, 2.0][(oid * 7 % 5) as usize];
            (oid, e_value, 30 + (oid * 13 % 11) as i32)
        })
        .collect();
    let mut hit_list = BlastnHitList::new(7);
    for &(oid, e_value, score) in &keys {
        hit_list.update(make_list(oid, &[(e_value, score)]));
    }
    hit_list.sort_by_evalue();
    let mut all: Vec<BlastnHspList> = keys
        .iter()
        .map(|&(oid, e_value, score)| {
            let mut list = make_list(oid, &[(e_value, score)]);
            list.update_best_evalue();
            list
        })
        .collect();
    all.sort_by(compare_hsp_lists);
    let best: Vec<u32> = all.iter().take(7).map(|list| list.oid).collect();
    assert_eq!(list_oids(&hit_list), best);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3267-3283
#[test]
fn hit_list_heapify_sorts_each_list_by_evalue() {
    let mut hit_list = BlastnHitList::new(1);
    hit_list.update(make_list(3, &[(1e-3, 10), (1e-20, 50)]));
    assert_eq!(hit_list.hsplist_array[0].hsps[0].raw_score, 10);
    hit_list.update(make_list(4, &[(1e-2, 70)]));
    assert!(hit_list.heapified);
    assert_eq!(list_oids(&hit_list), vec![3]);
    assert_eq!(hit_list.hsplist_array[0].hsps[0].raw_score, 50);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3172-3187
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:877-892
#[test]
fn hit_list_sort_purges_empty_lists_and_prune_caps_the_size() {
    let mut hit_list = BlastnHitList::new(10);
    hit_list.update(make_list(1, &[]));
    hit_list.update(make_list(2, &[(1e-5, 50)]));
    hit_list.update(make_list(3, &[]));
    hit_list.update(make_list(4, &[(1e-6, 60)]));
    for oid in 5..9 {
        hit_list.update(make_list(oid, &[(1e-4 * f64::from(oid), 40)]));
    }
    hit_list.sort_by_evalue();
    assert_eq!(list_oids(&hit_list), vec![4, 2, 5, 6, 7, 8]);
    hit_list.prune_by_size(4);
    assert_eq!(list_oids(&hit_list), vec![4, 2, 5, 6]);
    assert_eq!(hit_list.hsplist_count, 4);
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1264-1284
#[test]
fn test_blastn_outfmt7_header_matches_ncbi_shape() {
    let hit_lists = vec![None];
    let query_ids = vec![Arc::<[u8]>::from(&b"query"[..])];
    let subject_ids = vec![Some(Arc::<[u8]>::from(&b"subject"[..]))];
    let query_titles = vec![Arc::<[u8]>::from(&b"query full description"[..])];
    let mut output = Vec::new();

    write_output_blastn_hitlists_to_writer(
        &hit_lists,
        &mut output,
        &query_ids,
        &subject_ids,
        &crate::report::db_subject::HitLoads::default(),
        BlastnOutputFormat::TabularWithComments,
        &query_titles,
        "User specified sequence set (Input: subject.fasta)",
        &[],
        Some(1),
        None,
        None,
    )
    .unwrap();

    assert_eq!(
        String::from_utf8(output).unwrap(),
        concat!(
            "# BLASTN (LOSAT ",
            env!("CARGO_PKG_VERSION"),
            ")\n# Query: query full description\n# Database: User specified sequence set (Input: subject.fasta)\n# 0 hits found\n# BLAST processed 1 queries\n"
        )
    );
    // A query of a batch that NCBI did not search has no count.
    let mut output = Vec::new();
    write_output_blastn_hitlists_to_writer(
        &hit_lists,
        &mut output,
        &query_ids,
        &subject_ids,
        &crate::report::db_subject::HitLoads::default(),
        BlastnOutputFormat::TabularWithComments,
        &query_titles,
        "User specified sequence set (Input: subject.fasta)",
        &[true],
        None,
        None,
        None,
    )
    .unwrap();
    let text = String::from_utf8(output).unwrap();
    assert!(!text.contains("hits found"));
    // Without NCBI's epilog (an error in a later query batch), no last line.
    assert!(!text.contains("BLAST processed"));
}
