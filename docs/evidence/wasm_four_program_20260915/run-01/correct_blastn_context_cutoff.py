# NCBI reference: c++/src/algo/blast/core/blast_parameters.c:925-946,370-374
# searchsp = query_info->contexts[context].eff_searchsp;
# new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
# Apply the source-proven correction separately from optimization candidates.
import argparse
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('crate',type=Path);a=p.parse_args()
f=a.crate/'src/algorithm/blastn/blast_engine/run.rs';s=f.read_text()
old='''use super::super::ncbi_cutoffs::{
    compute_blastn_cutoff_score, compute_blastn_ungapped_params_from_score_freq,
    compute_eff_lengths_subject_mode_blastn, cutoff_score_max_from_evalue,
    GAP_TRIGGER_BIT_SCORE_NUCL,
};'''
new='''// NCBI reference: c++/src/algo/blast/core/blast_parameters.c:342-344,370-374
// gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 + kbp->logK) / kbp->Lambda);
// new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
use super::super::ncbi_cutoffs::{
    compute_blastn_ungapped_params_from_score_freq, cutoff_score_for_ungapped_extension,
    cutoff_score_max_from_evalue, gap_trigger_raw_score, GAP_TRIGGER_BIT_SCORE_NUCL,
};'''
assert s.count(old)==1;s=s.replace(old,new)
start=s.index('        let subject_len = s_len_full as i64;')
end=s.index('            cutoff_scores.push(cutoff);',start)
s=s[:start]+'''        let mut cutoff_scores: Vec<i32> = Vec::with_capacity(queries.len());
        let mut hit_saving_cutoff_scores: Vec<i32> = Vec::with_capacity(queries.len());
        // NCBI reference: c++/src/app/blast/blast_app_util.cpp:206-211;
        // c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-181;
        // c++/src/algo/blast/core/blast_engine.c:1434-1445
        // db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, true));
        // if (dbscan_mode) { ... m_iTotalLength += (Int8) (*iter)->length; }
        // if (db_length == 0) { BLAST_OneSubjectUpdateParameters(...); }
        // The CLI subject set has a nonzero total length. Its context search
        // space is retained across subjects, just as for output statistics.
        for query_idx in 0..queries.len() {
            // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:925-946
            // searchsp = query_info->contexts[context].eff_searchsp;
            // BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
            // params->cutoffs[context].cutoff_score_max = new_cutoff;
            // Both strand contexts share the same sequence length/search space.
            let hit_saving_cutoff = cutoff_score_max_from_evalue(
                evalue_threshold,
                query_eff_searchsp[query_idx * 2],
                &params_gapped_for_closure,
            );
            // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:342-374
            // gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 + kbp->logK) / kbp->Lambda);
            // new_cutoff = gap_trigger;
            // new_cutoff *= (Int4)sbp->scale_factor;
            // new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
            let gap_trigger = gap_trigger_raw_score(
                GAP_TRIGGER_BIT_SCORE_NUCL,
                &params_ungapped_for_closure,
            );
            let cutoff = cutoff_score_for_ungapped_extension(
                gap_trigger,
                hit_saving_cutoff,
                1.0,
            );
'''+s[end:]
f.write_text(s)
