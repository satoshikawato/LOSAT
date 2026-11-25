// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1450-1458
#[test]
fn pinned_official_subject_aa_cutoff_boundaries() {
    let golden = include_str!("../../../tests/unit/blastx_even_gap_cutoffs_expected.tsv");
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../docs/evidence/losatx_stage_a/inputs/generated");
    let mut options = super::super::stage_c_tests::options();
    options.gapped = false;
    for line in golden.lines() {
        let f: Vec<_> = line.split('\t').collect();
        let gap = f[0];
        assert_eq!(f[1], "B_LINK_SETUP");
        let records =
            super::super::input::read_fasta(&root.join(format!("S08_gap{gap}.fna")), false, false)
                .unwrap();
        let batch = super::super::query_setup::prepare_queries(&records, &options).unwrap();
        let mut parameters = super::super::parameters::score_block(&batch);
        let length: i32 = f[2].parse().unwrap();
        let db: i64 = f[3].parse().unwrap();
        super::super::parameters::subject_parameters(
            &batch,
            &mut parameters,
            &options,
            super::super::parameters::SubjectNumbers {
                total: db,
                count: 1,
                gumbel_total: db,
                shortest: length.max(10),
                average: db as i32,
                request_count: 1,
            },
            None,
        )
        .unwrap();
        let link = super::link_cutoffs(&batch, &parameters, length, db, None);
        assert_eq!(
            [link.cutoff_small_gap, link.cutoff_big_gap],
            [f[4].parse().unwrap(), f[5].parse().unwrap()],
            "{line}"
        );
        assert_eq!(
            link.gap_prob.to_bits(),
            u64::from_str_radix(f[6], 16).unwrap(),
            "{line}"
        );
    }
}
