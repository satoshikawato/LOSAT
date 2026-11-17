//! HSP Storage and Management
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c
//!
//! The filters of `s_FilterBlastResults` that every program shares (`-max_hsps`,
//! `-qcov_hsp_perc`), `Blast_HSPTest` (NCBI's `Int4` counts, widened to `i64`), and the frames
//! that `Blast_HSPGetTargetTranslation` never translates, with their refusal (TBLASTN, TBLASTX).

/// Whether an HSP fails NCBI's percent identity or minimal hit length test (`Blast_HSPTest`,
/// `s_HSPTest`: the HSP is deleted), for every program. The counts are NCBI's `Int4`s widened to
/// `i64`, so that the programs that count in `usize` pass their lengths unchanged (every `Int4`
/// is the same `i64` and the same `double`). The alignment length is negative for an ungapped
/// HSP whose end lies before its start (BLASTX's one-hit word finder's HSPs of length below 0
/// with `-window_size 0`, `aa_ungapped.c:1054,1083`): `align_length < min_hit_length` holds
/// there for every `min_hit_length` of 0 or more. The identity test multiplies the length by
/// `percent_identity` in `double`, as NCBI does. Neither test needs a guard for its disabled
/// value: with `-perc_identity` 0 (NCBI's range is 0 to 100, `blast_args.cpp:188-192`) and the
/// `min_hit_length` of 0 that NCBI's applications always leave (`calloc`, no option), both are
/// false for every length of 0 or more.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1446
/// ```c
///    *options = (BlastHitSavingOptions*) calloc(1, sizeof(BlastHitSavingOptions));
/// ```
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:993-1001,1027-1032
/// ```c
/// static Boolean s_HSPTest(const BlastHSP* hsp,
///                          const BlastHitSavingOptions* hit_options,
///                          Int4 align_length)
/// {
/// 	   return ((hsp->num_ident * 100.0 <
/// 			   align_length * hit_options->percent_identity) ||
/// 			   align_length < hit_options->min_hit_length) ;
///
/// }
/// ...
/// Boolean Blast_HSPTest(BlastHSP* hsp,
/// 		 	 	 	  const BlastHitSavingOptions* hit_options,
/// 		 	 	 	  Int4 align_length)
/// {
/// 	return s_HSPTest(hsp, hit_options, align_length);
/// }
/// ```
pub fn hsp_test(
    num_ident: i64,
    align_length: i64,
    percent_identity: f64,
    min_hit_length: i64,
) -> bool {
    // The counts are lengths of sequences in memory: far below 2^53, exact in `f64`.
    num_ident as f64 * 100.0 < align_length as f64 * percent_identity
        || align_length < min_hit_length
}

/// The query coverage of an HSP in percent: its query extent (`offset` to `end`, end
/// exclusive, in the letters of its context) over the length of its context, plus 0.5 below
/// 99 (`Blast_HSPGetQueryCoverage`).
///
/// NCBI reference: c++/src/algo/blast/core/blast_hits.c:1034-1044
/// ```c
/// double Blast_HSPGetQueryCoverage(const BlastHSP* hsp, Int4 query_length)
/// {
/// 	double pct = 0;
///     if(query_length > 0) {
///             pct = 100.0 * (double) (hsp->query.end - hsp->query.offset)/ (double) query_length;
///             if(pct < 99)
///             	pct +=0.5;
///     }
///     return pct;
/// }
/// ```
pub fn hsp_query_coverage(offset: i64, end: i64, query_length: i64) -> f64 {
    let mut pct = 0.0;
    if query_length > 0 {
        pct = 100.0 * (end - offset) as f64 / query_length as f64;
        if pct < 99.0 {
            pct += 0.5;
        }
    }
    pct
}

/// Whether an HSP covers less than `min_pct` percent of its query context
/// (`Blast_HSPQueryCoverageTest`: the HSP is deleted).
///
/// NCBI reference: c++/src/algo/blast/core/blast_hits.c:1046-1051
/// ```c
/// Boolean Blast_HSPQueryCoverageTest(BlastHSP* hsp,
///         						   double min_query_coverage_pct,
///         						   Int4 query_length)
/// {
///      double hsp_coverage = Blast_HSPGetQueryCoverage( hsp, query_length);
///      return (hsp_coverage < min_query_coverage_pct);
/// }
/// ```
pub fn hsp_query_coverage_test(min_pct: f64, offset: i64, end: i64, query_length: i64) -> bool {
    hsp_query_coverage(offset, end, query_length) < min_pct
}

/// `Blast_HSPListReapByQueryCoverage` over the HSPs of one list, in their order: the HSPs
/// whose query coverage (`hsp_query_coverage_test` of the `extent` of each: query offset,
/// query end and the length of its context) is below `min_pct` are removed. A `min_pct` of
/// 0 is NCBI's unset option and removes nothing. Returns whether an HSP was removed (the
/// caller then sets the list's best e-value again, and frees a list left empty).
///
/// NCBI reference: c++/src/algo/blast/core/blast_hits.c:2010-2045
/// ```c
///    if ((hsp_list == NULL) || (hsp_list->hspcnt == 0) || (hit_options->query_cov_hsp_perc == 0))
///       return 0;
///
///    hsp_array = hsp_list->hsp_array;
///    for (index = 0; index < hsp_list->hspcnt; index++) {
///       hsp = hsp_array[index];
///       ASSERT(hsp != NULL);
///       if ( Blast_HSPQueryCoverageTest(hsp, hit_options->query_cov_hsp_perc,
/// 			                          query_info->contexts[hsp->context].query_length)) {
///          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
///          removed = TRUE;
///       } else {
///          if (index > hsp_cnt)
///             hsp_array[hsp_cnt] = hsp_array[index];
///          hsp_cnt++;
///       }
///    }
///
///    hsp_list->hspcnt = hsp_cnt;
///
///    if (removed) {
///        hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
///    }
/// ```
pub fn reap_by_query_coverage<T>(
    hsps: &mut Vec<T>,
    min_pct: f64,
    extent: impl Fn(&T) -> (i64, i64, i64),
) -> bool {
    if hsps.is_empty() || min_pct == 0.0 {
        return false;
    }
    let count = hsps.len();
    hsps.retain(|hsp| {
        let (offset, end, query_length) = extent(hsp);
        !hsp_query_coverage_test(min_pct, offset, end, query_length)
    });
    hsps.len() != count
}

/// `Blast_TrimHSPListByMaxHsps`: the first `max_hsps` HSPs of a list, in its order, are
/// kept (0 keeps all). The list's best e-value is not set again.
///
/// NCBI reference: c++/src/algo/blast/core/blast_hits.c:2049-2069
/// ```c
/// Int2 Blast_TrimHSPListByMaxHsps(BlastHSPList* hsp_list,
///                                 const BlastHitSavingOptions* hit_options)
/// {
///    ...
///    if ((hsp_list == NULL) ||
/// 	   (hit_options->max_hsps_per_subject == 0) ||
/// 	   (hsp_list->hspcnt <= hit_options->max_hsps_per_subject))
///       return 0;
///
///    hsp_max = hit_options->max_hsps_per_subject;
///    hsp_array = hsp_list->hsp_array;
///    for (index = hsp_max; index < hsp_list->hspcnt; index++) {
///       hsp_array[index] = Blast_HSPFree(hsp_array[index]);
///    }
///
///    hsp_list->hspcnt = hsp_max;
///    return 0;
/// }
/// ```
pub fn trim_by_max_hsps<T>(hsps: &mut Vec<T>, max_hsps: usize) {
    if max_hsps == 0 || hsps.len() <= max_hsps {
        return;
    }
    hsps.truncate(max_hsps);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1284;
// c++/src/algo/blast/core/blast_hits.c:1161-1163,1233-1234
// ```c
//       retval->translations = (Uint1**) calloc(num_frames, sizeof(Uint1*));
// ...
//     /* skip translation if full translation has already been done */
//     if (target_t->partial && (start ||
//         (stop < target_t->subject_blk->length / CODON_LENGTH -3)))
// ...
//     /* +1 as the first byte is a sentinel. */
//     return target_t->translations[context] - target_t->range[2*context] + 1;
// ```
/// Whether NCBI's `Blast_HSPGetTargetTranslation` never translates a frame of a subject of
/// `subject_length` bases (the searched range, `-subject_loc`: NCBI's `subject_blk->length`): its
/// test is false at the first call (`start` and `stop` 0) for fewer than 12 bases and stays so,
/// and it returns the translation it never made (`calloc`'d to NULL) plus 1, from which every
/// caller reads next. TBLASTN's traceback reads (`tblastn::search_gapped::TargetTranslation`) and
/// TBLASTX's ungapped re-evaluation (`blast_hits.c:2700-2705`) refuse there
/// (`x_untranslated_frame_refusal`).
#[inline]
pub(crate) fn target_translation_never_made(subject_length: usize) -> bool {
    subject_length / 3 <= 3
}

/// The refusal of `program`'s `step` (the NCBI step that reads the translation) of an alignment
/// on a subject whose frames NCBI never translates (`target_translation_never_made`): NCBI BLAST+
/// 2.17.0 reads from address 1 on (NULL + 1) and dies with SIGSEGV (PD-LOSAT-NCBI-DEFECTS 1.46 for
/// TBLASTN, 1.49 for TBLASTX, rule 3).
#[cold]
#[inline(never)]
pub(crate) fn x_untranslated_frame_refusal(program: &str, step: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "{program}: {step} of an alignment on a subject (or its -subject_loc range) of fewer than 12 bases, where NCBI BLAST+ reads the translation of the subject's frame, which it never makes for a subject this short, and crashes; LOSAT does not support it"
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn query_coverage_adds_a_half_below_99() {
        // 95 of 1000 letters: 9.5 % + 0.5.
        assert_eq!(hsp_query_coverage(5, 100, 1000), 10.0);
        // 988 of 1000: 98.8 % + 0.5 = 99.3; 992 of 1000 stays 99.2.
        assert_eq!(hsp_query_coverage(0, 988, 1000), 98.8 + 0.5);
        assert_eq!(hsp_query_coverage(0, 992, 1000), 99.2);
        assert_eq!(hsp_query_coverage(0, 10, 0), 0.0);
    }

    #[test]
    fn reap_keeps_the_order_and_reports_removals() {
        let mut hsps = vec![(0, 50), (0, 95), (10, 100), (0, 9)];
        let extent = |&(offset, end): &(i64, i64)| (offset, end, 1000);
        assert!(!reap_by_query_coverage(&mut hsps, 0.0, extent));
        assert_eq!(hsps.len(), 4);
        // Coverage 5.5, 10.0, 9.5, 1.4: a limit of 9.5 keeps the second and the third.
        assert!(reap_by_query_coverage(&mut hsps, 9.5, extent));
        assert_eq!(hsps, vec![(0, 95), (10, 100)]);
        assert!(!reap_by_query_coverage(&mut hsps, 9.5, extent));
    }

    #[test]
    fn hsp_test_deletes_a_negative_length_at_the_default_options() {
        // BLASTX's options: percent_identity 0 and min_hit_length 0. A negative length (an
        // ungapped HSP whose end lies before its start) is deleted by the length test (its
        // identity test is `0.0 < -0.0`, false); a length of 0 is kept: both tests are false
        // (`0.0 < 0.0`, `0 < 0`).
        for length in [-4, -1] {
            assert!(hsp_test(0, length, 0.0, 0));
        }
        assert!(!hsp_test(0, 0, 0.0, 0));
        assert!(!hsp_test(3, 7, 0.0, 0));
        // `-perc_identity -0` (in NCBI's range): `0.0 < -0.0` and `500.0 < -0.0` are false.
        assert!(!hsp_test(0, 5, -0.0, 0));
        assert!(!hsp_test(5, 5, -0.0, 0));
        // The identity test: 3 of 7 is 42.9 %, below 50 and not below 40.
        assert!(hsp_test(3, 7, 50.0, 0));
        assert!(!hsp_test(3, 7, 40.0, 0));
        // The minimal hit length.
        assert!(hsp_test(7, 7, 0.0, 8));
        assert!(!hsp_test(7, 7, 0.0, 7));
    }

    /// NCBI's C API (`Blast_HSPTest` called directly on the pinned NCBI build) on the
    /// boundaries: identity equal to the limit kept and just below deleted, the length test,
    /// and the empty HSP at the default options.
    #[test]
    fn hsp_test_matches_the_ncbi_api_oracle() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/hsp_test_20260924/oracle.tsv"
        );
        let oracle = std::fs::read_to_string(path).unwrap();
        for line in oracle.lines() {
            let f: Vec<_> = line.split('\t').collect();
            let observed = hsp_test(
                f[1].parse().unwrap(),
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
            );
            assert_eq!(observed, f[5] == "1", "{line}");
        }
        assert_eq!(oracle.lines().count(), 6);
    }

    /// No NCBI counterpart: NCBI's `Int4` cannot hold these counts. The `usize` counts of
    /// BLASTN, TBLASTN and TBLASTX are NCBI's `Int4` widened to `i64`: a length past `Int4` is
    /// tested as itself, not capped.
    #[test]
    fn hsp_test_takes_lengths_past_int4() {
        let long = i64::from(i32::MAX) + 10;
        assert!(!hsp_test(long, long, 100.0, 0));
        assert!(hsp_test(long - 1, long, 100.0, 0));
        assert!(!hsp_test(0, long, 0.0, long));
        assert!(hsp_test(0, long - 1, 0.0, long));
    }

    /// NCBI's test `stop < length/3 - 3` with `stop` 0 (an `Int4` division, negative below 9
    /// bases) holds from 12 bases on.
    #[test]
    fn a_target_translation_is_never_made_below_12_bases() {
        for length in 0..60usize {
            let ncbi_translates = 0 < length as i32 / 3 - 3;
            assert_eq!(
                target_translation_never_made(length),
                !ncbi_translates,
                "{length}"
            );
        }
        assert!(target_translation_never_made(11));
        assert!(!target_translation_never_made(12));
        assert_eq!(
            x_untranslated_frame_refusal("TBLASTX", "the ungapped re-evaluation").to_string(),
            "TBLASTX: the ungapped re-evaluation of an alignment on a subject (or its -subject_loc range) of fewer than 12 bases, where NCBI BLAST+ reads the translation of the subject's frame, which it never makes for a subject this short, and crashes; LOSAT does not support it"
        );
    }

    #[test]
    fn trim_keeps_the_first_hsps() {
        let mut hsps = vec![3, 1, 2];
        trim_by_max_hsps(&mut hsps, 0);
        assert_eq!(hsps, vec![3, 1, 2]);
        trim_by_max_hsps(&mut hsps, 3);
        assert_eq!(hsps, vec![3, 1, 2]);
        trim_by_max_hsps(&mut hsps, 2);
        assert_eq!(hsps, vec![3, 1]);
    }
}
