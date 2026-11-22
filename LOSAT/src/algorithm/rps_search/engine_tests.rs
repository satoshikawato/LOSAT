use super::*;

// C11 6.3.1.4: a double converts to `int` truncated toward zero when the truncated value is
// representable; otherwise (and for NaN) the behaviour is undefined, which `int4_of` refuses.
#[test]
fn int4_conversions_truncate_toward_zero_inside_32_bits_only() {
    assert_eq!(int4_of(2_147_483_647.9), Some(i32::MAX));
    assert_eq!(int4_of(2_147_483_648.0), None);
    assert_eq!(int4_of(-2_147_483_648.9), Some(i32::MIN));
    assert_eq!(int4_of(-2_147_483_649.0), None);
    assert_eq!(int4_of(-0.9), Some(0));
    assert_eq!(int4_of(-2.596), Some(-2));
    assert_eq!(int4_of(64.9), Some(64));
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e300, -1e300] {
        assert_eq!(int4_of(value), None, "{value}");
    }
}

// blast_parameters.c:216-221: `scale_factor` times the ceiling of `x_dropoff` times `LN2` over
// `Lambda`, converted to `Int4`: the ceiling before the scale, one conversion of the product.
#[test]
fn the_ungapped_x_drop_is_the_ceiling_times_the_scale() {
    let lambda = 0.3176;
    // 7 bits: ceil(15.28) = 16 raw, 1,600 at scale 100.
    assert_eq!(ungapped_x_dropoff(7.0, lambda, 1.0), Some(16));
    assert_eq!(ungapped_x_dropoff(7.0, lambda, 100.0), Some(1600));
    // A tiny X-drop is still one raw unit.
    assert_eq!(ungapped_x_dropoff(1e-300, lambda, 200.0), Some(200));
    // ceil(x) * scale just inside and just outside 32 bits.
    let bits_for = |raw: f64| raw * lambda / std::f64::consts::LN_2;
    assert_eq!(
        ungapped_x_dropoff(bits_for(21_474_836.0) - 1e-6, lambda, 100.0),
        Some(2_147_483_600)
    );
    assert_eq!(
        ungapped_x_dropoff(bits_for(21_474_837.0), lambda, 100.0),
        None
    );
    assert_eq!(ungapped_x_dropoff(1e300, lambda, 1.0), None);
}

/// An `RpsRun` of no query with the given X-drops in bits, BLOSUM62 11/1's gapped block and
/// no composition-based statistics.
fn x_drop_run<'a>(
    gumbel: &'a BlastGumbelBlk,
    options: &'a RpsSearch,
    gapped_bits: f64,
    final_bits: f64,
) -> RpsRun<'a> {
    RpsRun {
        queries: Vec::new(),
        gapped: KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        },
        gumbel: Some(gumbel),
        cutoff_subject_length: 0,
        evalue: 10.0,
        prelim_evalue: 10.0,
        cbs: 0,
        sbp: None,
        hitlist_size: 500,
        prelim_hitlist_size: 500,
        window: 40,
        // rps_mini's rows: 190 + 180 + 190 letters and 3 sentinel rows.
        concat_length: 563,
        sorts_results: true,
        prelim_spaces: None,
        x_drop_ungapped_bits: 7.0,
        x_drop_gapped_bits: gapped_bits,
        x_drop_final_bits: final_bits,
        gap_trigger_bits: 22.0,
        options,
        program: "RPSBLAST",
    }
}

/// A valid query of `len` letters (between two sentinels) with BLOSUM62's ungapped block.
fn cutoff_query(letters: &[u8], search_space: i64) -> RpsQuery<'_> {
    let kbp = KarlinParams {
        lambda: 0.3176,
        k: 0.134,
        h: 0.4012,
        alpha: 0.0,
        beta: 0.0,
    };
    RpsQuery {
        frames: vec![RpsFrame {
            frame: 0,
            aa_seq: letters,
            aa_seq_nomask: None,
            valid: true,
            kbp_std: kbp,
            search_space,
            first_query_length: letters.len() - 2,
        }],
        batch_kbp: kbp,
    }
}

// blast_parameters.c:922-945 and 978-990 (cited at `RpsCutoffs`): with the Gumbel block the
// hit cutoff is `BLAST_SpougeEtoS` of the preliminary e-value, without it
// `BLAST_Cutoffs` of `kbp_gap` and `eff_searchsp`, both times the scale (100 here).
#[test]
fn the_hit_cutoff_reads_the_gumbel_block_or_the_search_space_at_the_scale() {
    use crate::config::ProteinScoringSpec;
    let letters = vec![1u8; 302];
    let options = RpsSearch::default();
    let gumbel =
        crate::stats::lookup_protein_gumbel_params(&ProteinScoringSpec::default(), 3_234_567)
            .unwrap();
    let mut run = x_drop_run(&gumbel, &options, 15.0, 25.0);
    run.cutoff_subject_length = 250;
    let query = cutoff_query(&letters, 1_000_000);
    let with_block = run
        .cutoffs(&query.frames[0], query.frames[0].search_space, 1, 100.0)
        .unwrap();
    assert_eq!(
        with_block.hit_cutoff,
        crate::stats::blast_spouge_etos(10.0, &run.gapped, &gumbel, 300, 250) * 100
    );
    // Without the block: ceil(ln(K * searchsp / E) / Lambda) = ceil(ln(4100) / 0.267) = 32.
    run.gumbel = None;
    let without = run
        .cutoffs(&query.frames[0], query.frames[0].search_space, 1, 100.0)
        .unwrap();
    assert_eq!(without.hit_cutoff, 3_200);
    // The search space is the query's, not a constant.
    let other = cutoff_query(&letters, 1_000_000_000);
    assert_eq!(
        run.cutoffs(&other.frames[0], other.frames[0].search_space, 1, 100.0)
            .unwrap()
            .hit_cutoff,
        5_800
    );
    // The word cutoff is the gap trigger, at most the hit cutoff (blast_parameters.c:370).
    assert!(without.word_cutoff <= without.hit_cutoff);
    // The preliminary e-value is `-evalue` as given (the stretch is the caller's).
    run.prelim_evalue = 1e300;
    assert_eq!(
        run.cutoffs(&query.frames[0], query.frames[0].search_space, 1, 100.0)
            .unwrap()
            .hit_cutoff,
        100
    );
}

// blast_parameters.c:455-469 (cited at `RpsRun::gapped_x_drops`): each X-drop truncated
// toward zero, the final one the larger of its own and the preliminary one (compared as
// doubles), both times the scale above 1; the values NCBI's probes give (DESIGN_ps11c §2).
#[test]
fn gapped_x_drops_are_ncbis_conversions_at_the_scale() {
    let gumbel = BlastGumbelBlk {
        lambda: 0.267,
        c: 0.0,
        g: 0.0,
        a: 0.0,
        alpha: 0.0,
        sigma: 0.0,
        a_un: 0.0,
        alpha_un: 0.0,
        b: 0.0,
        beta: 0.0,
        tau: 0.0,
        db_length: 0,
    };
    // The X-drops as given on the command line (the refusals name them).
    let given = |gapped: f64, last: f64| RpsSearch {
        x_drop_gapped_bits: Some(gapped),
        x_drop_final_bits: Some(last),
        ..RpsSearch::default()
    };
    // BLOSUM62 11/1 at the scale.
    let drops = |gapped: f64, last: f64, scale: f64| {
        let options = given(gapped, last);
        x_drop_run(&gumbel, &options, gapped, last)
            .gapped_x_drops(scale, 11 * scale as i32, scale as i32)
            .map_err(|error| error.to_string())
    };
    // The defaults: 15 bits = 38 raw, 25 bits = 64 raw.
    assert_eq!(drops(15.0, 25.0, 1.0), Ok((38, 64)));
    assert_eq!(drops(15.0, 25.0, 100.0), Ok((3800, 6400)));
    // 0 and negative gapped X-drops are kept as they are (the kernels raise them).
    assert_eq!(drops(0.0, 25.0, 1.0), Ok((0, 64)));
    assert_eq!(drops(-1.0, -1.0, 100.0), Ok((-200, -200)));
    assert_eq!(drops(-1.0, 25.0, 100.0), Ok((-200, 6400)));
    // The final X-drop is at least the preliminary one: 40 and 10 bits give 103 and 103.
    assert_eq!(drops(40.0, 10.0, 1.0), Ok((103, 103)));
    // A final X-drop far below 32 bits takes the preliminary one: no undefined conversion.
    assert_eq!(drops(15.0, -1e300, 100.0), Ok((3800, 3800)));
    for (gapped, last, scale, option) in [
        (1e300, 25.0, 1.0, "-xdrop_gap of 1e300 bits"),
        (-1e300, 25.0, 1.0, "-xdrop_gap of -1e300 bits"),
        (15.0, 1e300, 1.0, "-xdrop_gap_final of 1e300 bits"),
        // 1e7 bits = 25,960,568 raw: inside 32 bits at scale 1, outside at scale 100.
        (1e7, 25.0, 100.0, "-xdrop_gap of 10000000.0 bits"),
        (15.0, 1e7, 100.0, "-xdrop_gap_final of 10000000.0 bits"),
    ] {
        let error = drops(gapped, last, scale).unwrap_err();
        assert!(
            error.contains(option) && error.contains(&format!("scale {scale}")),
            "{error}"
        );
    }
    assert_eq!(drops(1e7, 1e7, 1.0), Ok((25_960_568, 25_960_568)));
    // The kernels' window: 4.1e8 bits = 1,064,383,310 raw, below 2^30 - 103 over a gap
    // extension of 1; 4.2e8 bits = 1,090,343,879 raw, above it: refused for the preliminary
    // X-drop, and for the final one where the CBS 0 traceback runs (this run, no redo). At
    // scale 2 the product with the scale already leaves 32 bits.
    assert!(drops(4.1e8, 25.0, 1.0).is_ok());
    assert!(drops(15.0, 4.1e8, 1.0).is_ok());
    for scale in [1.0, 2.0] {
        let error = drops(4.2e8, 25.0, scale).unwrap_err();
        assert!(error.contains("-xdrop_gap of 420000000.0 bits"), "{error}");
        let error = drops(15.0, 4.2e8, scale).unwrap_err();
        assert!(
            error.contains("-xdrop_gap_final of 420000000.0 bits"),
            "{error}"
        );
    }
    // NCBI's defaults never leave 32 bits but at an extreme scale, which the refusal names.
    let defaults = RpsSearch::default();
    let error = x_drop_run(&gumbel, &defaults, 15.0, 25.0)
        .gapped_x_drops(1e8, 1_100_000_000, 100_000_000)
        .unwrap_err()
        .to_string();
    assert!(
        error.starts_with("a profile database of the scale 100000000 (its scaled gap costs"),
        "{error}"
    );
}

// blast_gapalign.c:426-429, 457-464, 641-644 (cited at `kernel_x_drop_fits`): the X-drop
// raised to the gap costs, over the gap extension, plus 103 is the DP window; twice a window
// of 2^30 leaves `Int4`.
#[test]
fn the_kernels_window_of_an_x_drop_stays_in_32_bits() {
    assert_eq!(
        crate::core::x_drop::KERNEL_X_DROP_QUOTIENT_LIMIT,
        1_073_741_721
    );
    assert!(kernel_x_drop_fits(1_073_741_720, 11, 1));
    assert!(!kernel_x_drop_fits(1_073_741_721, 11, 1));
    // Scale 2: the quotient truncates.
    assert!(kernel_x_drop_fits(2_147_483_441, 22, 2));
    assert!(!kernel_x_drop_fits(2_147_483_442, 22, 2));
    // From scale 3 on no `Int4` X-drop reaches it.
    assert!(kernel_x_drop_fits(i32::MAX, 33, 3));
    // An X-drop below the gap costs is raised to them; negative ones too.
    assert!(kernel_x_drop_fits(0, 11, 1));
    assert!(kernel_x_drop_fits(i32::MIN, 11, 1));
}

// blast_traceback.c:1241-1261: with composition-based statistics only the redo runs after
// the preliminary stage, whose gap extension is 32 times the database's: the final X-drop
// then never reaches the kernels' window limit, the preliminary one still does.
#[test]
fn the_final_x_drop_reaches_the_kernels_only_without_the_redo() {
    static FREQ_RATIOS: [[f64; crate::core::blast_util::BLASTAA_SIZE];
        crate::core::blast_util::BLASTAA_SIZE] =
        [[0.0; crate::core::blast_util::BLASTAA_SIZE]; crate::core::blast_util::BLASTAA_SIZE];
    static BACKGROUND: [f64; 20] = [0.05; 20];
    let sbp = crate::core::psi::PsiScoreBlk {
        matrix: [[0; crate::core::blast_util::BLASTAA_SIZE]; crate::core::blast_util::BLASTAA_SIZE],
        ideal_lambda: 0.3176,
        ideal_k: 0.134,
        gap_std_k: 0.041,
        gap_psi: crate::core::psi::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
        },
        freq_ratios: &FREQ_RATIOS,
        bit_scale_factor: 2,
        background: &BACKGROUND,
    };
    let gumbel = BlastGumbelBlk {
        lambda: 0.267,
        c: 0.0,
        g: 0.0,
        a: 0.0,
        alpha: 0.0,
        sigma: 0.0,
        a_un: 0.0,
        alpha_un: 0.0,
        b: 0.0,
        beta: 0.0,
        tau: 0.0,
        db_length: 0,
    };
    let options = RpsSearch {
        x_drop_gapped_bits: Some(15.0),
        x_drop_final_bits: Some(4.2e8),
        ..RpsSearch::default()
    };
    let mut run = x_drop_run(&gumbel, &options, 15.0, 4.2e8);
    run.cbs = 1;
    run.sbp = Some(&sbp);
    // A final X-drop of 4.2e8 bits leaves 32 bits at the redo's Lambda: the redo's check.
    let error = run.gapped_x_drops(1.0, 11, 1).unwrap_err().to_string();
    assert!(error.contains("an -xdrop_gap_final of"), "{error}");
    // (At the redo's Lambda a final X-drop is about 32 times its raw value, so one above the
    // kernels' window limit is always refused by this conversion.) 1.5e7 bits are
    // 1,246,107,290 at the redo's Lambda: accepted.
    run.x_drop_final_bits = 1.5e7;
    assert!(run.gapped_x_drops(1.0, 11, 1).is_ok());
    // The preliminary X-drop is checked at every mode.
    run.x_drop_gapped_bits = 4.2e8;
    let error = run.gapped_x_drops(1.0, 11, 1).unwrap_err().to_string();
    assert!(error.contains("an -xdrop_gap of"), "{error}");
}

// blast_kappa.c:2381-2383 (cited at `RedoRun::new`): the redo's X-drop at Lambda over the
// scaling factor 32; 25 bits give 2,076 raw, below the scaled final X-drop 6,400 of scale 100.
#[test]
fn the_redo_x_drop_is_the_larger_of_the_bits_at_the_redo_lambda_and_the_final_one() {
    let lambda = 0.267 / SCALING_FACTOR;
    assert_eq!(
        int4_of(RedoRun::x_dropoff_of(25.0, lambda, 6400)),
        Some(6400)
    );
    assert_eq!(int4_of(RedoRun::x_dropoff_of(25.0, lambda, 64)), Some(2076));
    assert_eq!(
        int4_of(RedoRun::x_dropoff_of(-1.0, lambda, -200)),
        Some(-83)
    );
    // 2.5e7 bits = 2,076,845,484 raw: inside; 1e8 bits: outside.
    assert!(int4_of(RedoRun::x_dropoff_of(2.5e7, lambda, 64)).is_some());
    assert_eq!(int4_of(RedoRun::x_dropoff_of(1e8, lambda, 64)), None);
    assert_eq!(int4_of(RedoRun::x_dropoff_of(-1e300, lambda, 64)), Some(64));
}

#[test]
fn the_first_query_of_a_batch_gets_the_k_of_the_highest_mark() {
    let marks = [Some(2), None, Some(5), None, None, Some(1)];
    let karlin_k = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06];
    let k = batch_gapped_k(&[0..3, 3..5, 5..6], &marks, &karlin_k, 0.041, false);
    assert_eq!(k[0], RPS_K_MULT * 0.06);
    assert_eq!(k[1], 0.041);
    assert_eq!(k[2], 0.041);
    assert_eq!(k[3], 0.041);
    assert_eq!(k[5], RPS_K_MULT * 0.02);
    // Read by score (rpsblast with composition-based statistics): the mark of the last query
    // that has one, whatever its ordinal id.
    let marks = [Some(5), None, Some(2), None, None, Some(1)];
    let k = batch_gapped_k(&[0..3, 3..5, 5..6], &marks, &karlin_k, 0.041, true);
    assert_eq!(k[0], RPS_K_MULT * 0.03);
    assert_eq!(k[3], 0.041);
    assert_eq!(k[5], RPS_K_MULT * 0.02);
    let k = batch_gapped_k(&[0..3, 3..5, 5..6], &marks, &karlin_k, 0.041, false);
    assert_eq!(k[0], RPS_K_MULT * 0.06);
}

fn hsp(
    oid: u32,
    e_value: f64,
    q_start: usize,
    q_end: usize,
) -> crate::algorithm::blastp::hsp::BlastpHsp {
    let mut hsp = crate::algorithm::blastp::hsp::BlastpHsp::from_hit(crate::common::Hit {
        identity: 0.0,
        length: q_end - q_start + 1,
        mismatch: 0,
        gapopen: 0,
        q_start,
        q_end,
        s_start: 1,
        s_end: 10,
        e_value,
        bit_score: 10.0,
        num_ident: 0,
        query_frame: 0,
        query_length: 100,
        q_idx: 0,
        s_idx: oid,
        raw_score: 20,
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: None,
        num_positives: 0,
    });
    hsp.e_value = e_value;
    hsp
}

// tabular.cpp:973-1007 and showalign.cpp:2120-2135 (cited at `formatter_identities`): the
// reports count the identities on the query as read and the profile's consensus, across the
// gaps of the edit script (after `swap_roles`: Del steps the profile, Ins the query).
#[test]
fn the_reports_identities_count_the_letters_as_read() {
    use crate::algorithm::blastp::hsp::{BlastpHitList, BlastpHspList};
    use crate::common::GapEditOp;
    // Profile 1 starts at 6 of the concatenated consensus (profile 0 holds 5 letters + 1).
    let enc = |text: &str| -> Vec<u8> {
        text.bytes()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect()
    };
    let mut consensus = enc("AAAAA");
    consensus.push(0);
    consensus.extend(enc("MKTAYWWIAKQ"));
    let offsets = [0, 6, 18];
    // The query as read (its searched letters could have had X over KTAY, a lower-case run;
    // the traceback's count of 5 below is that of the masked letters).
    let letters = enc("GMKTAYIAKQ");
    // The letters as read between two sentinels (`RpsFrame::aa_seq_nomask`).
    let mut nomask = vec![0u8];
    nomask.extend(&letters);
    nomask.push(0);
    let query = RpsQuery {
        frames: vec![RpsFrame {
            frame: 0,
            aa_seq: &nomask,
            aa_seq_nomask: Some(&nomask),
            valid: true,
            kbp_std: KarlinParams::default(),
            search_space: 0,
            first_query_length: letters.len(),
        }],
        batch_kbp: KarlinParams::default(),
    };
    // Query 2..10 against profile 1..11: MKTAY / MKTAY, a gap of WW in the query, IAKQ / IAKQ.
    let mut aligned = hsp(1, 1e-5, 2, 10);
    aligned.s_start = 1;
    aligned.s_end = 11;
    aligned.length = 11;
    aligned.num_ident = 5;
    aligned.mismatch = 4;
    aligned.gap_info = Some(vec![
        GapEditOp::Sub(5),
        GapEditOp::Del(2),
        GapEditOp::Sub(4),
    ]);
    let mut hit_list = BlastpHitList::new(10);
    hit_list.update(BlastpHspList {
        oid: 1,
        query_index: 0,
        hsps: vec![aligned],
        best_evalue: 1e-5,
    });
    formatter_identities(&mut hit_list, &query, &consensus, &offsets);
    let counted = &hit_list.hsplist_array[0].hsps[0];
    assert_eq!(counted.num_ident, 9);
    assert_eq!(counted.mismatch, 0);
    assert_eq!(counted.identity, 900.0 / 11.0);
    // An Ins steps the query alone.
    let mut inserted = hsp(1, 1e-5, 1, 6);
    inserted.s_start = 1;
    inserted.s_end = 5;
    inserted.length = 6;
    inserted.gap_info = Some(vec![GapEditOp::Ins(1), GapEditOp::Sub(5)]);
    let mut hit_list = BlastpHitList::new(10);
    hit_list.update(BlastpHspList {
        oid: 1,
        query_index: 0,
        hsps: vec![inserted],
        best_evalue: 1e-5,
    });
    formatter_identities(&mut hit_list, &query, &consensus, &offsets);
    assert_eq!(hit_list.hsplist_array[0].hsps[0].num_ident, 5);
}

#[test]
fn a_coverage_filter_that_removes_a_best_hsp_reorders_and_cuts_the_volumes_list() {
    // The batch's first query has 100 letters.
    let letters = [0u8; 102];
    let query = RpsQuery {
        frames: vec![RpsFrame {
            frame: 0,
            aa_seq: &letters,
            aa_seq_nomask: None,
            valid: true,
            kbp_std: KarlinParams::default(),
            search_space: 0,
            first_query_length: 100,
        }],
        batch_kbp: KarlinParams::default(),
    };
    use crate::algorithm::blastp::hsp::{BlastpHitList, BlastpHspList};
    let list = |oid: u32, hsps: Vec<_>, best: f64| BlastpHspList {
        oid,
        query_index: 0,
        hsps,
        best_evalue: best,
    };
    let build = || {
        let mut hit_list = BlastpHitList::new(10);
        // Subject 1: its best HSP covers 10 % of the 100-letter query, its second 80 %.
        hit_list.update(list(
            1,
            vec![hsp(1, 1e-10, 1, 10), hsp(1, 1e-3, 11, 90)],
            1e-10,
        ));
        hit_list.update(list(2, vec![hsp(2, 1e-5, 1, 90)], 1e-5));
        hit_list
    };
    let oids = |hit_list: &BlastpHitList| -> Vec<u32> {
        hit_list.hsplist_array[..hit_list.hsplist_count]
            .iter()
            .map(|list| list.oid)
            .collect()
    };
    let mut hit_list = build();
    hit_list.sort_by_evalue();
    assert_eq!(oids(&hit_list), vec![1, 2]);
    // 50 % coverage drops subject 1's best HSP: its best E-value becomes 1e-3, after 2's.
    let mut hit_list = build();
    finish_volume_list(&mut hit_list, 0, 50.0, &query, 10, true);
    assert_eq!(oids(&hit_list), vec![2, 1]);
    // And the list is cut to the hit list size after the sort.
    let mut hit_list = build();
    finish_volume_list(&mut hit_list, 0, 50.0, &query, 1, true);
    assert_eq!(oids(&hit_list), vec![2]);
    // Without a filter the order is the list's own.
    let mut hit_list = build();
    finish_volume_list(&mut hit_list, 0, 0.0, &query, 10, true);
    assert_eq!(oids(&hit_list), vec![1, 2]);
    // A unit whose total length is 0 or negative keeps the order of the filtered list
    // (blast_traceback.c:1768-1771): subject 1 stays first, then the cut.
    let mut hit_list = build();
    finish_volume_list(&mut hit_list, 0, 50.0, &query, 10, false);
    assert_eq!(oids(&hit_list), vec![1, 2]);
    let mut hit_list = build();
    finish_volume_list(&mut hit_list, 0, 50.0, &query, 1, false);
    assert_eq!(oids(&hit_list), vec![1]);
}

/// The files of the fixture database `rps_mini` (tests/db_search/rps/, three profiles)
/// inflated into a directory of this test, and the database opened.
fn rps_mini(name: &str) -> (std::path::PathBuf, crate::seqdb::SeqDb) {
    rps_mini_without(name, &[])
}

/// `rps_mini` without the files of the extensions `left_out`.
fn rps_mini_without(name: &str, left_out: &[&str]) -> (std::path::PathBuf, crate::seqdb::SeqDb) {
    use std::io::Read;
    let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/rps");
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-rps-engine-{name}-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    for extension in [
        "aux", "loo", "rps", "pin", "psq", "phr", "pdb", "pos", "pot", "ptf", "pto",
    ] {
        if left_out.contains(&extension) {
            continue;
        }
        let mut bytes = Vec::new();
        flate2::read::GzDecoder::new(
            std::fs::read(fixtures.join(format!("rps_mini.{extension}.gz")))
                .unwrap()
                .as_slice(),
        )
        .read_to_end(&mut bytes)
        .unwrap();
        std::fs::write(dir.join(format!("rps_mini.{extension}")), bytes).unwrap();
    }
    let db = crate::seqdb::SeqDb::open(
        dir.join("rps_mini").to_str().unwrap(),
        crate::seqdb::SeqType::Protein,
    )
    .unwrap();
    (dir, db)
}

/// NCBIstdaa letters of `text` between two sentinels.
fn encoded(text: &str) -> Vec<u8> {
    const LETTERS: &[u8] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
    let mut letters = vec![0u8];
    letters.extend(text.bytes().map(|letter| {
        LETTERS
            .iter()
            .position(|&l| l == letter.to_ascii_uppercase())
            .unwrap_or(21) as u8
    }));
    letters.push(0);
    letters
}

// D2 of PS12s: a volume's search with every valid context of every query a unit of the pool
// (the query finished by the worker that searches its last context) gives each query the
// results of the search one context after another on one worker; a query whose cutoffs are
// refused gets that refusal as its result, the others theirs.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn units_of_contexts_give_the_results_of_the_search_in_order() {
    use crate::algorithm::rps_search::volumes::{Consensus, RpsVolumeSet, RpsVolumeSetup};
    let (dir, db) = rps_mini("units");
    let volumes = RpsVolumeSet::of(&db);
    let volume = &volumes.volumes[0];
    let mut setup = RpsVolumeSetup::new("RPSTBLASTN");
    assert!(setup.check(volume, false).unwrap().is_ok());
    let Ok(files) = setup.reopen(volume, false).unwrap() else {
        panic!("rps_mini opens");
    };
    // The profiles' fragments (tests/db_search/rps/build.sh) as the frames of queries: windows
    // of the WSSV proteins with the profiles' letters, and two queries of other letters.
    let proteins = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/rps/q40.faa"),
    )
    .unwrap();
    let sequences: Vec<String> = proteins
        .split('>')
        .filter(|record| !record.is_empty())
        .map(|record| record.lines().skip(1).collect())
        .collect();
    let windows = [
        (32, 0, 80),
        (32, 60, 190),
        (18, 300, 400),
        (18, 380, 480),
        (4, 1100, 1200),
        (4, 1180, 1290),
        (0, 0, 120),
        (1, 10, 90),
    ];
    let frames: Vec<Vec<u8>> = (0..18)
        .map(|i| {
            let (record, start, end) = windows[i % windows.len()];
            let sequence = &sequences[record];
            encoded(&sequence[start.min(sequence.len())..end.min(sequence.len())])
        })
        .collect();
    let kbp = KarlinParams {
        lambda: 0.3176,
        k: 0.134,
        h: 0.4012,
        alpha: 0.0,
        beta: 0.0,
    };
    let options = RpsSearch::default();
    // BLOSUM62 11/1's Gumbel block (the database's 1,000,000 letters).
    let gumbel = crate::stats::spouge::lookup_protein_gumbel_params(
        &crate::config::ProteinScoringSpec {
            matrix: crate::config::ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
        },
        1_000_000,
    )
    .unwrap();
    // Six queries of six frames (one of them invalid), the second with a Lambda that makes its
    // ungapped X-drop leave 32 bits (refused: NCBI's default X-drop, so the refusal names the
    // database's scale).
    let queries: Vec<RpsQuery<'_>> = (0..6)
        .map(|q| RpsQuery {
            frames: (0..NUM_FRAMES)
                .map(|f| {
                    let aa_seq = &frames[(q * NUM_FRAMES + f) % frames.len()];
                    RpsFrame {
                        frame: [1, 2, 3, -1, -2, -3][f],
                        aa_seq,
                        aa_seq_nomask: None,
                        valid: f != 4,
                        kbp_std: KarlinParams {
                            lambda: if q == 1 { 1e-300 } else { kbp.lambda },
                            ..kbp
                        },
                        search_space: 20_000_000,
                        first_query_length: aa_seq.len() - 2,
                    }
                })
                .collect(),
            batch_kbp: kbp,
        })
        .collect();
    let run = RpsRun {
        queries,
        gapped: KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        },
        gumbel: Some(&gumbel),
        cutoff_subject_length: 180,
        evalue: 10.0,
        prelim_evalue: 10.0,
        cbs: 0,
        sbp: None,
        hitlist_size: 500,
        prelim_hitlist_size: 500,
        window: 40,
        // rps_mini's rows: 190 + 180 + 190 letters and 3 sentinel rows.
        concat_length: 563,
        sorts_results: true,
        prelim_spaces: None,
        x_drop_ungapped_bits: 7.0,
        x_drop_gapped_bits: 15.0,
        x_drop_final_bits: 25.0,
        gap_trigger_bits: 22.0,
        options: &options,
        program: "RPSTBLASTN",
    };
    let consensus = Consensus::new(&setup, volume, &files);
    let shown = |results: RpsVolumeResults| -> Vec<String> {
        assert!(results.late.is_none());
        results
            .refused
            .into_iter()
            .zip(results.searched)
            .map(|(refusal, result)| match (refusal, result) {
                (Some(refusal), None) => format!("refused: {refusal}"),
                (None, Some(Ok(traced))) => {
                    format!("{:?} {:?}", traced.hit_list, traced.k_mark)
                }
                (None, Some(Err(error))) => format!("failed: {error}"),
                _ => unreachable!("a query is refused or searched"),
            })
            .collect()
    };
    let serial = shown(
        run.search_volume(
            &files,
            &consensus,
            &crate::utils::threading::SearchPool::x_serial(),
            &RpsWorkers::new(1),
        )
        .unwrap(),
    );
    assert_eq!(serial.len(), 6);
    assert!(
        serial[1].starts_with("refused: a profile database of the scale 100"),
        "{}",
        serial[1]
    );
    let with_hits = serial
        .iter()
        .filter(|shown| shown.starts_with("Some(BlastpHitList"))
        .count();
    assert!(with_hits >= 4, "{serial:?}");
    for threads in [2, 4] {
        let workers = RpsWorkers::new(threads);
        let parallel = crate::utils::threading::with_search_pool(threads, "RPSTBLASTN", |pool| {
            // Twice on the same workers: the kept scratch gives the same results.
            let first = shown(run.search_volume(&files, &consensus, pool, &workers)?);
            let again = shown(run.search_volume(&files, &consensus, pool, &workers)?);
            assert_eq!(first, again);
            Ok(first)
        })
        .unwrap();
        assert_eq!(parallel, serial, "{threads} threads");
    }
    let _ = std::fs::remove_dir_all(dir);
}

// S1 and N1 of the review of PS12s: a batch's refusal of a query's cutoffs comes before the
// search errors of the batch's earlier queries, and the volume's refusal of its Lambda at the
// scale after the cutoffs (NCBI sets the cutoffs of the whole batch before its scan,
// blast_engine.c:1372-1374, cited at `search_volume`; its e-values come after). Query 0 has
// hits, so its traceback reads the consensus letters, whose file (.psq) is missing: its search
// fails; query 1's ungapped X-drop leaves 32 bits: refused at its cutoffs.
#[test]
fn a_batchs_cutoff_refusal_comes_before_its_search_errors() {
    use crate::algorithm::rps_search::volumes::{Consensus, RpsVolumeSet, RpsVolumeSetup};
    let proteins = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/rps/q40.faa"),
    )
    .unwrap();
    let record: String = proteins
        .split('>')
        .filter(|record| !record.is_empty())
        .nth(32)
        .unwrap()
        .lines()
        .skip(1)
        .collect();
    let letters = encoded(&record[..80]);
    let options = RpsSearch::default();
    let gumbel = crate::stats::spouge::lookup_protein_gumbel_params(
        &crate::config::ProteinScoringSpec::default(),
        1_000_000,
    )
    .unwrap();
    let queries = || {
        let mut queries = vec![
            cutoff_query(&letters, 20_000_000),
            cutoff_query(&letters, 20_000_000),
        ];
        queries[1].frames[0].kbp_std.lambda = 1e-300;
        queries
    };
    let run = |lambda: f64| RpsRun {
        queries: queries(),
        gapped: KarlinParams {
            lambda,
            ..x_drop_run(&gumbel, &options, 15.0, 25.0).gapped
        },
        cutoff_subject_length: 180,
        ..x_drop_run(&gumbel, &options, 15.0, 25.0)
    };
    let search = |left_out: &[&str], lambda: f64| {
        let (dir, db) = rps_mini_without("cutoffs", left_out);
        let volumes = RpsVolumeSet::of(&db);
        let volume = &volumes.volumes[0];
        let mut setup = RpsVolumeSetup::new("RPSBLAST");
        assert!(setup.check(volume, false).unwrap().is_ok());
        let Ok(files) = setup.reopen(volume, false).unwrap() else {
            panic!("rps_mini opens");
        };
        let consensus = Consensus::new(&setup, volume, &files);
        let results = run(lambda)
            .search_volume(
                &files,
                &consensus,
                &crate::utils::threading::SearchPool::x_serial(),
                &RpsWorkers::new(1),
            )
            .unwrap();
        let _ = std::fs::remove_dir_all(dir);
        results
    };
    let errors = |results: RpsVolumeResults, spans: &[usize]| -> Vec<String> {
        results
            .batches(spans)
            .into_iter()
            .map(|batch| match batch {
                Ok(traced) => format!("{} lists", traced.len()),
                Err(error) => error.to_string(),
            })
            .collect()
    };
    let lambda = 0.267;
    // With every file: query 0 has a hit list, query 1 is refused.
    let whole = search(&[], lambda);
    assert!(whole.refused[1].is_some() && whole.refused[0].is_none());
    assert!(matches!(&whole.searched[0], Some(Ok(traced)) if traced.hit_list.is_some()));
    let scale_refusal = "a profile database of the scale 100 (its scaled";
    assert!(errors(whole, &[2])[0].starts_with(scale_refusal));
    // Without the consensus letters: query 0's search fails; in one batch query 1's refusal
    // still comes first, in a batch of its own query 0's error ends the batches.
    let missing = |spans: &[usize]| errors(search(&["psq"], lambda), spans);
    assert!(
        missing(&[2])[0].starts_with(scale_refusal),
        "{:?}",
        missing(&[2])
    );
    let alone = missing(&[1, 1]);
    assert_eq!(alone.len(), 1);
    assert!(alone[0].contains("rps_mini.psq"), "{alone:?}");
    // A Lambda that the scale does not give back (0.4389734947748931 / 100 * 100 differs):
    // the volume's refusal comes after the batch's cutoff refusals, before any search.
    let late = 0.438_973_494_774_893_1;
    assert_ne!(late / 100.0 * 100.0, late);
    let lambda_refusal = "a profile database of the scale 100 (NCBI BLAST+ divides Lambda";
    let refused = search(&[], late);
    assert!(refused.searched.iter().all(Option::is_none));
    assert!(errors(refused, &[2])[0].starts_with(scale_refusal));
    let first = errors(search(&[], late), &[1, 1]);
    assert_eq!(first.len(), 1);
    assert!(first[0].starts_with(lambda_refusal), "{first:?}");
    // `all` reads the run's queries as one batch (rpsblast's).
    let Err(error) = search(&[], lambda).all() else {
        panic!("the run's queries are refused");
    };
    assert!(error.to_string().starts_with(scale_refusal), "{error}");
}
