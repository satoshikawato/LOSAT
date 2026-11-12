//! The RPS-BLAST traceback without composition-based statistics (`s_RPSComputeTraceback` with
//! `compositionBasedStats` 0): for each preliminary list of a query, the profile's rows rescaled
//! to the query's composition (`RPSRescalePssm`), the profile's K, `Blast_TracebackFromHSPList`
//! with the profile as the query and the user query as the subject (no containment test before
//! each alignment), then `s_HSPListPostTracebackUpdate`: the roles swapped back, the e-values at
//! the database's scale, the scores divided back and the bit scores.

use anyhow::Result;

use super::engine::RpsQuery;
use super::redo::RpsQueryStats;
use super::rescale::{residue_probabilities, rps_rescale_pssm};
use super::ProfileSlice;
use super::RPS_K_MULT;
use crate::algorithm::blastp::hsp::{BlastpHitList, BlastpHsp, BlastpHspList};
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::pssm_search::gapalign::PssmGapScratch;
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::algorithm::pssm_search::traceback0::{trace_hsps, TraceHsps, TraceSubject};
use crate::common::GapEditOp;
use crate::stats::BlastGumbelBlk;
use crate::stats::KarlinParams;

/// What the traceback of every query of a volume reads besides the query.
pub(crate) struct RpsTraceback<'a> {
    /// The concatenated rows and the profiles' offsets.
    pub(crate) rows: &'a [[i32; 28]],
    pub(crate) offsets: &'a [i32],
    /// The consensus letters of the profiles, at the profiles' offsets (the sentinel rows hold
    /// any byte).
    pub(crate) consensus: &'a [u8],
    /// The database's alphabet (26 or 28 columns).
    pub(crate) alphabet: usize,
    /// The `.aux` scale (`score_params->scale_factor`).
    pub(crate) scale: f64,
    /// The `.aux` K of each profile.
    pub(crate) karlin_k: &'a [f64],
    /// The gap costs and the final X-drop at the scale.
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    pub(crate) x_drop_final: i32,
    /// `kbp_gap_std` (BLOSUM62 11/1's table values), whose K each list replaces.
    pub(crate) gapped: KarlinParams,
    /// `hit_options->expect_value`.
    pub(crate) expect_value: f64,
    /// `hit_options->hitlist_size`.
    pub(crate) hitlist_size: usize,
    /// The composition-based redo (modes 1 to 3), else the rescaled PSSM's traceback.
    pub(crate) redo: Option<super::redo::RpsRedoParams<'a>>,
}

/// The memory of one worker's traceback.
pub(crate) struct RpsTraceScratch {
    pub(super) gap: PssmGapScratch,
    pub(super) pssm: QueryPssm,
    /// The frequency ratios of the profile being redone.
    pub(super) ratios: Vec<[f64; 28]>,
    /// The bytes of those rows as the `.freq` file holds them.
    pub(super) freq_bytes: Vec<u8>,
}

impl RpsTraceScratch {
    pub(crate) fn new() -> Self {
        Self {
            gap: PssmGapScratch::new(),
            pssm: QueryPssm { rows: Vec::new() },
            ratios: Vec::new(),
            freq_bytes: Vec::new(),
        }
    }
}

/// The final hit list of a query and the highest profile whose rescaled PSSM it wrote the K of
/// (the K the query's batch may print, `engine.rs`).
pub(crate) struct RpsTraced {
    pub(crate) hit_list: Option<BlastpHitList>,
    pub(crate) k_mark: Option<u32>,
}

/// `s_BlastHSPListRPSUpdate`'s part before the sort: the query and the subject of each HSP
/// switched back (the HSPs are in the profile's and the user query's local coordinates; the
/// frames too: a translated query's frame becomes the HSP's query frame), the insertions and
/// deletions of its edit script switched, and the user query's fields (`query_length` of the
/// HSP's query frame: the length of the user query or of its frame).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:165-184
/// ```c
///    hsp = hsplist->hsp_array;
///    for (i = 0; i < hsplist->hspcnt; i++) {
///
///       /* switch query and subject (which are already in local coordinates) */
///       tmp = hsp[i]->query;
///       hsp[i]->query = hsp[i]->subject;
///       hsp[i]->subject = tmp;
///
///       /* Change the traceback information to reflect the query and subject
///          sequences getting switched */
///       s_BlastHSPRPSUpdate(hsp[i]);
///
///
///       /* If query was nucleotide, set context, because it is needed in order
///          to pick correct Karlin block for calculating bit scores. There are
///          6 contexts corresponding to each nucleotide query sequence. */
///       if (program == eBlastTypeRpsTblastn) {
///           hsp[i]->context = BLAST_FrameToContext(hsp[i]->query.frame, program);
///       }
///    }
///    Blast_HSPListSortByScore(hsplist);
/// ```
/// (The sort is the caller's.)
pub(super) fn swap_roles(hsp_list: &mut BlastpHspList, query_length: impl Fn(i32) -> usize) {
    for hsp in &mut hsp_list.hsps {
        std::mem::swap(&mut hsp.q_start, &mut hsp.s_start);
        std::mem::swap(&mut hsp.q_end, &mut hsp.s_end);
        std::mem::swap(&mut hsp.gapped_q_start, &mut hsp.gapped_s_start);
        std::mem::swap(&mut hsp.query_frame, &mut hsp.subject_frame);
        if let Some(ops) = hsp.gap_info.as_mut() {
            for op in ops.iter_mut() {
                *op = match *op {
                    GapEditOp::Ins(n) => GapEditOp::Del(n),
                    GapEditOp::Del(n) => GapEditOp::Ins(n),
                    sub => sub,
                };
            }
        }
        hsp.query_length = query_length(hsp.query_frame);
        hsp.query_context = hsp_list.query_index as i32;
        hsp.q_idx = hsp_list.query_index;
        hsp.s_idx = hsp_list.oid;
    }
}

/// The identities of a volume's final HSPs of one query as NCBI's formatter counts them: the
/// columns where the query's letters as read (`query_letters`, without the SEG and lower-case
/// masks and without the X of `Blast_MaskUnsupportedAA`) equal the profile's consensus letters.
/// The traceback counts on the searched letters (`num_ident` of `final_hsp` and of the redo, as
/// NCBI's core does), which differ where the query is masked; NCBI's reports do not read that
/// count but the aligned sequences of their scope (the outfmt 0 `Identities`, the tabular
/// `nident`, `pident` and `mismatch`), so the HSPs of a run whose writer prints the HSP's count
/// (outfmt 0, the default outfmt 6 and 7) take these; a tabular report with a field list counts on
/// its rows and keeps the traceback's count, which opens NCBI's alignment block for `nident`
/// (`blastp::tabular::needs_rows`, `rps_hooks::rps_search`). The roles are swapped already
/// (`swap_roles`: query offsets on the user query, subject offsets on the profile).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:973-1007
/// ```c++
///     if (x_IsFieldRequested(eQuerySeq) ||
///         ...
///         (x_IsFieldRequested(eNumIdentical) && !kNoFetchSequence) ||
///         (x_IsFieldRequested(eMismatches) && !kNoFetchSequence) ||
///         (x_IsFieldRequested(ePercentIdentical) && !kNoFetchSequence)) {
///
///         alnVec->SetGapChar('-');
///         ...
///         alnVec->GetWholeAlnSeqString(0, m_QuerySeq);
///         alnVec->GetWholeAlnSeqString(1, m_SubjectSeq);
///         ...
///             num_ident = 0;
///             ...
///             for (unsigned int i = 0;
///                  i < min(m_QuerySeq.size(), m_SubjectSeq.size());
///                  ++i) {
///                 if (m_QuerySeq[i] == m_SubjectSeq[i]) {
///                     ++num_ident;
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2120-2135
/// ```c++
/// void CDisplaySeqalign::x_FillIdentityInfo(const string& sequence_standard,
///                                           const string& sequence ,
///                                           int& match, int& positive,
///                                           string& middle_line)
/// {
///     match = 0;
///     ...
///     for(int i=0; i<min_length; i++){
///         if(sequence_standard[i]==sequence[i]){
///             ...
///             match ++;
/// ```
pub(super) fn formatter_identities(
    hit_list: &mut BlastpHitList,
    query: &RpsQuery<'_>,
    consensus: &[u8],
    offsets: &[i32],
) {
    for hsp_list in hit_list
        .hsplist_array
        .iter_mut()
        .take(hit_list.hsplist_count)
    {
        let profile = &consensus[offsets[hsp_list.oid as usize] as usize..];
        for hsp in &mut hsp_list.hsps {
            let Some(ops) = hsp.gap_info.as_ref() else {
                continue;
            };
            // The letters as read of the HSP's query context (its frame).
            let Some(nomask) = query.context(hsp.query_frame).aa_seq_nomask else {
                continue;
            };
            let query_letters = &nomask[1..nomask.len() - 1];
            let (mut q, mut s) = (hsp.q_start - 1, hsp.s_start - 1);
            let mut num_ident = 0usize;
            for op in ops {
                match *op {
                    GapEditOp::Sub(n) => {
                        let n = n as usize;
                        num_ident += query_letters[q..q + n]
                            .iter()
                            .zip(&profile[s..s + n])
                            .filter(|(a, b)| a == b)
                            .count();
                        q += n;
                        s += n;
                    }
                    GapEditOp::Del(n) => s += n as usize,
                    GapEditOp::Ins(n) => q += n as usize,
                }
            }
            hsp.mismatch += hsp.num_ident as i64 - num_ident as i64;
            hsp.num_ident = num_ident;
            hsp.identity = if hsp.length > 0 {
                (num_ident as f64 * 100.0) / (hsp.length as f64)
            } else {
                0.0
            };
        }
    }
}

impl RpsTraceback<'_> {
    /// The traceback of the query of letters `query` (without sentinels) and index `q_idx` from
    /// its preliminary hit list.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1191-1236
    /// ```c
    ///          else {
    ///
    ///              /* for no CBS option rescale PSSM matrix with the half-way CBS */
    ///              /* replace the PSSM and the Karlin values for this DB sequence
    ///                 and this query sequence. */
    ///              sbp->psi_matrix->pssm->data =
    ///                  RPSRescalePssm(score_params->scale_factor,
    ///                                 one_query->length, one_query->sequence,
    ///                                 seq_arg.seq->length,
    ///                                 rpsblast_pssms + db_seq_start, sbp);
    ///              /* The composition of the query could have caused this one
    ///                 subject sequence to produce a bad PSSM. This should
    ///                 not be a fatal error, so just go on to the next subject
    ///                 sequence */
    ///              if (sbp->psi_matrix->pssm->data == NULL) {
    ///                  ...
    ///                  hsp_list = Blast_HSPListFree(hsp_list);
    ///                  BlastSeqSrcReleaseSequence(seq_src, &seq_arg);
    ///                  continue;
    ///              }
    ///          }
    ///
    ///          sbp->kbp_gap[0]->K = RPS_K_MULT * karlin_k[hsp_list->oid];
    ///          sbp->kbp_gap[0]->logK = log(RPS_K_MULT * karlin_k[hsp_list->oid]);
    ///       }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1272-1305
    /// ```c
    ///       if (hsp_list->hspcnt == 0) {
    ///          hsp_list = Blast_HSPListFree(hsp_list);
    ///          continue;
    ///       }
    ///
    ///       /* Save this HSP list in the results structure. */
    ///       Blast_HSPResultsInsertHSPList(results, hsp_list,
    ///                                     hit_params->options->hitlist_size);
    ///    }
    ///    ...
    ///    Blast_HSPResultsSortByEvalue(results);
    /// ```
    /// The lists are read by ordinal id (the stream of a batch without composition-based
    /// statistics sorts them so); a list's traceback reads nothing another writes, and the hit
    /// list's order is total, so the read order does not change the hit list.
    pub(crate) fn trace_query(
        &self,
        prelim: BlastpHitList<XPrelimHsp>,
        query: &RpsQuery<'_>,
        q_idx: u32,
        gumbel: Option<&BlastGumbelBlk>,
        scratch: &mut RpsTraceScratch,
    ) -> Result<RpsTraced> {
        let first = &query.frames[0];
        // A translated query's contexts are its frames; a protein query has one, frame 0.
        let translated = first.frame != 0;
        let letters = &first.aa_seq[1..first.aa_seq.len().max(2) - 1];
        if let Some(redo) = self.redo.as_ref() {
            let statistics = RpsQueryStats {
                batch_kbp: &query.batch_kbp,
                stats: first.stats(gumbel),
            };
            return self.trace_query_redo(redo, prelim, letters, q_idx, statistics, scratch);
        }
        let probabilities = residue_probabilities(letters);
        // The frames' letters without their sentinels, in the context order (NCBI's
        // `translation_buffer + frame_offsets[context] + 1`).
        let frames: Vec<&[u8]> = query
            .frames
            .iter()
            .map(|frame| &frame.aa_seq[1.min(frame.aa_seq.len())..frame.aa_seq.len().max(1) - 1])
            .collect();
        let mut lists: Vec<BlastpHspList<XPrelimHsp>> = prelim
            .hsplist_array
            .into_iter()
            .take(prelim.hsplist_count)
            .collect();
        lists.sort_by_key(|list| list.oid);
        let mut hit_list: Option<BlastpHitList> = None;
        let mut k_mark: Option<u32> = None;
        for list in lists {
            let oid = list.oid as usize;
            let start = self.offsets[oid] as usize;
            let length = (self.offsets[oid + 1] - self.offsets[oid] - 1) as usize;
            let consensus = &self.consensus[start..start + length];
            let hsps = list.hsps.iter().map(|hsp| hsp.hsp);
            let (hsps, karlin, stat_length) = if translated {
                // The profile's rows at the database's scale, without the rescale and the K of
                // composition-based statistics 0 (blast_traceback.c:1184-1198); each HSP aligns
                // on its frame after the containment test (blast_traceback.c:402-417).
                let rows = ProfileSlice {
                    rows: self.rows,
                    start,
                    len: length,
                };
                let (hsps, stat_length) = trace_hsps(
                    &TraceHsps {
                        rows: &rows,
                        query_nomask: consensus,
                        query_length: length,
                        tree_query_bound: length as i32 + 1,
                        subject: TraceSubject::Frames {
                            frames: &frames,
                            // NCBI reference (598d8ae6):
                            // c++/src/algo/blast/core/blast_traceback.c:371-373 (the subject's
                            // length is 0 for a translated subject: `subject_blk->length`, the
                            // frames with their sentinels, over `CODON_LENGTH`)
                            // ```c
                            //    tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
                            //                                  0, (subject_length > 0 ? subject_length :
                            //                                  subject_blk->length / CODON_LENGTH) + 1);
                            // ```
                            tree_bound: (query.translation_length() / 3) as i32 + 1,
                        },
                        gap_open: self.gap_open,
                        gap_extend: self.gap_extend,
                        x_drop_final: self.x_drop_final,
                        contained_before: true,
                    },
                    hsps,
                    &mut scratch.gap,
                );
                (hsps, self.gapped, stat_length)
            } else {
                let profile_rows = &self.rows[start..start + length];
                if !rps_rescale_pssm(
                    self.scale,
                    &probabilities,
                    profile_rows,
                    self.alphabet,
                    &mut scratch.pssm.rows,
                )? {
                    continue;
                }
                k_mark = Some(list.oid);
                let (hsps, _) = trace_hsps(
                    &TraceHsps {
                        rows: &scratch.pssm,
                        query_nomask: consensus,
                        query_length: length,
                        tree_query_bound: length as i32 + 1,
                        subject: TraceSubject::One(letters),
                        gap_open: self.gap_open,
                        gap_extend: self.gap_extend,
                        x_drop_final: self.x_drop_final,
                        contained_before: false,
                    },
                    hsps,
                    &mut scratch.gap,
                );
                let karlin = KarlinParams {
                    k: RPS_K_MULT * self.karlin_k[oid],
                    ..self.gapped
                };
                (hsps, karlin, length)
            };
            let mut hsp_list = BlastpHspList {
                oid: list.oid,
                query_index: q_idx,
                hsps,
                best_evalue: f64::MAX,
            };
            self.post_update(&mut hsp_list, &karlin, query, gumbel, length, stat_length);
            if hsp_list.hsps.is_empty() {
                continue;
            }
            hit_list
                .get_or_insert_with(|| BlastpHitList::new(self.hitlist_size))
                .update(hsp_list);
        }
        if let Some(hit_list) = hit_list.as_mut() {
            hit_list.sort_by_evalue();
        }
        Ok(RpsTraced { hit_list, k_mark })
    }

    /// `s_HSPListPostTracebackUpdate` of one list: the roles switched back, then each HSP's
    /// e-value (`m` the user query's length, `n` the profile's; for a translated query `m` the
    /// length of the frame of the list's last HSP aligned, the traceback's `stat_length`, and `n`
    /// the length of the HSP's frame, whose search space it reads without the Gumbel block), the
    /// reap by `-evalue`, the scores rescaled and the bit scores.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1876-1883 (quoted at
    /// `EvalueStats::e_value`: `isRPS` for both RPS programs, `subject_length` the traceback's
    /// `stat_length`, `query_info` the user query's, whose contexts are a translated query's
    /// frames)
    fn post_update(
        &self,
        hsp_list: &mut BlastpHspList,
        karlin: &KarlinParams,
        query: &RpsQuery<'_>,
        gumbel: Option<&BlastGumbelBlk>,
        profile_length: usize,
        stat_length: usize,
    ) {
        let translated = query.frames[0].frame != 0;
        swap_roles(hsp_list, |frame| query.context(frame).length());
        sort_by_score(&mut hsp_list.hsps);
        let scaled = KarlinParams {
            lambda: karlin.lambda / self.scale,
            ..*karlin
        };
        for hsp in &mut hsp_list.hsps {
            let context = query.context(hsp.query_frame);
            hsp.e_value = if translated {
                context.stats(gumbel).e_value(
                    hsp.raw_score,
                    &scaled,
                    stat_length as i32,
                    context.length() as i32,
                )
            } else {
                context.stats(gumbel).e_value(
                    hsp.raw_score,
                    &scaled,
                    context.length() as i32,
                    profile_length as i32,
                )
            };
        }
        hsp_list.hsps.retain(|hsp| hsp.e_value <= self.expect_value);
        for hsp in &mut hsp_list.hsps {
            hsp.raw_score = ((f64::from(hsp.raw_score) + 0.5 * self.scale) / self.scale) as i32;
        }
        sort_by_score(&mut hsp_list.hsps);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1926-1928
        // ```c
        //       hsp->bit_score =
        //          (hsp->score*kbp[hsp->context]->Lambda -
        //           kbp[hsp->context]->logK) / NCBIMATH_LN2;
        // ```
        let log_k = karlin.k.ln();
        for hsp in &mut hsp_list.hsps {
            hsp.bit_score =
                (f64::from(hsp.raw_score) * karlin.lambda - log_k) / std::f64::consts::LN_2;
        }
        hsp_list.best_evalue = hsp_list
            .hsps
            .iter()
            .map(|hsp| hsp.e_value)
            .fold(f64::MAX, f64::min);
    }
}

/// `Blast_HSPListSortByScore` on final HSPs: sorted with `ScoreCompareHSPs` (score down, subject
/// offset up, subject end down, query offset up, query end down) when not sorted already, stable
/// as glibc's qsort.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1350
/// ```c
///     if (0 == (result = BLAST_CMP(h2->score, h1->score)) &&
///         0 == (result = BLAST_CMP(h1->subject.offset, h2->subject.offset)) &&
///         0 == (result = BLAST_CMP(h2->subject.end, h1->subject.end)) &&
///         0 == (result = BLAST_CMP(h1->query.offset, h2->query.offset))) {
///         result = BLAST_CMP(h2->query.end, h1->query.end);
///     }
/// ```
pub(super) fn sort_by_score(hsps: &mut [BlastpHsp]) {
    let compare = |a: &BlastpHsp, b: &BlastpHsp| {
        b.raw_score
            .cmp(&a.raw_score)
            .then_with(|| a.s_start.cmp(&b.s_start))
            .then_with(|| b.s_end.cmp(&a.s_end))
            .then_with(|| a.q_start.cmp(&b.q_start))
            .then_with(|| b.q_end.cmp(&a.q_end))
    };
    if hsps
        .windows(2)
        .any(|pair| compare(&pair[0], &pair[1]) == std::cmp::Ordering::Greater)
    {
        hsps.sort_by(compare);
    }
}

/// `s_FilterBlastResults` for RPS-BLAST: `-max_hsps` (the first HSPs of each list) and
/// `-qcov_hsp_perc`, whose query length is that of the context of the HSP, 0 for every RPS-BLAST
/// HSP: the length of the first query of the batch (`first_query_length`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:845-866
/// ```c
///       for (subject_index = 0;
///            subject_index < hit_list->hsplist_count; ++subject_index) {
///     	  	  BlastHSPList * hsp_list = hit_list->hsplist_array[subject_index];
///     	  	  if(hit_options->max_hsps_per_subject) {
///     	  		  Blast_TrimHSPListByMaxHsps(hsp_list, hit_options);
///     	  	  }
///     	  	  if(hit_options->query_cov_hsp_perc) {
///     	  		  Blast_HSPListReapByQueryCoverage(hsp_list, hit_options, query_info, program_number);
///     	  		  if(hsp_list->hspcnt == 0){
///     	  			hit_list->hsplist_array[subject_index] =  Blast_HSPListFree(hsp_list);
///     	  		  }
///     	  	  }
///          ...
///        }
///        if(hit_options->query_cov_hsp_perc) {
///             Blast_HitListPurgeNullHSPLists(hit_list);
///        }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2026-2044
/// ```c
///    for (index = 0; index < hsp_list->hspcnt; index++) {
///       hsp = hsp_array[index];
///       ASSERT(hsp != NULL);
///       if ( Blast_HSPQueryCoverageTest(hsp, hit_options->query_cov_hsp_perc,
/// 			                          query_info->contexts[hsp->context].query_length)) {
///          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
///          removed = TRUE;
///       ...
///    if (removed) {
///        hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
///    }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1034-1051
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
///
/// Boolean Blast_HSPQueryCoverageTest(BlastHSP* hsp,
///         						   double min_query_coverage_pct,
///         						   Int4 query_length)
/// {
///      double hsp_coverage = Blast_HSPGetQueryCoverage( hsp, query_length);
///      return (hsp_coverage < min_query_coverage_pct);
/// }
/// ```
pub(crate) fn filter_results(
    hit_list: &mut BlastpHitList,
    max_hsps: usize,
    query_cov_hsp_perc: f64,
    first_query_length: impl Fn(i32) -> usize,
) {
    if max_hsps == 0 && query_cov_hsp_perc <= 0.0 {
        return;
    }
    for hsp_list in hit_list
        .hsplist_array
        .iter_mut()
        .take(hit_list.hsplist_count)
    {
        if max_hsps > 0 && hsp_list.hsps.len() > max_hsps {
            hsp_list.hsps.truncate(max_hsps);
        }
        if query_cov_hsp_perc > 0.0 {
            let count = hsp_list.hsps.len();
            hsp_list.hsps.retain(|hsp| {
                let first_query_length = first_query_length(hsp.query_frame);
                let mut pct = 0.0;
                if first_query_length > 0 {
                    pct = 100.0 * (hsp.q_end + 1 - hsp.q_start) as f64 / first_query_length as f64;
                    if pct < 99.0 {
                        pct += 0.5;
                    }
                }
                pct >= query_cov_hsp_perc
            });
            if hsp_list.hsps.len() != count {
                hsp_list.best_evalue = hsp_list
                    .hsps
                    .iter()
                    .map(|hsp| hsp.e_value)
                    .fold(f64::MAX, f64::min);
            }
        }
    }
    if query_cov_hsp_perc > 0.0 {
        hit_list.purge_empty_hsp_lists();
    }
}
