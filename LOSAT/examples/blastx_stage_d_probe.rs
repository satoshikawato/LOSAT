//! Comparison-only diagnostic; no public BLASTX search or report success.
use anyhow::{bail, Result};
use LOSAT::{
    algorithm::blastx::{
        input::read_fasta,
        linking::LinkedHspList,
        search::search_preliminary,
        statistics::{evaluate_gapped, reap_preliminary},
    },
    cli::{Cli, Commands},
};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1802-1810
// ```c++
//     /* Sort the HSP array by score */
//     Blast_HSPListSortByScore(hsp_list);
//
//     /* Find and fill the best e-value */
//     hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
//     for (index = 1; index < hsp_list->hspcnt; ++index) {
//         if (hsp_list->hsp_array[index]->evalue < hsp_list->best_evalue)
//             hsp_list->best_evalue = hsp_list->hsp_array[index]->evalue;
//     }
// ```
fn dump(stage: &str, call: usize, oid: usize, list: &LinkedHspList) {
    if list.hsps.is_empty() {
        return;
    }
    println!("{stage}_COUNT\t{call}\t{oid}\t{}", list.hsps.len());
    for (i, h) in list.hsps.iter().enumerate() {
        let p = &h.hsp;
        println!(
            "{stage}\t{call}\t{oid}\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:016x}",
            p.context,
            p.score,
            p.frame,
            p.q_start,
            p.q_end,
            p.q_gapped_start,
            p.s_start,
            p.s_end,
            p.s_gapped_start,
            0,
            h.num,
            h.evalue.to_bits()
        );
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-899
// ```c++
//     if (hit_params->link_hsp_params) {
//         status = BLAST_LinkHsps(program_number, hsp_list_out, query_info,
//                   subject->length, gap_align->sbp, hit_params->link_hsp_params,
//                   score_options->gapped_calculation);
//     } else if (!Blast_ProgramIsPhiBlast(program_number)
//            && !(isRPS && !sbp->gbp)
//            /* do not calculate E-values for mapping */
//            && program_number != eBlastTypeMapping ) {
//         /* Calculate e-values for all HSPs. Skip this step
//            for PHI or RPS with old FSC, since calculating the E values
//            requires precomputation that has not been done yet */
//         double scale_factor = 1.0;
//         if (isRPS) {
//             scale_factor = score_params->scale_factor;
//         }
//         Blast_HSPListGetEvalues(program_number, query_info,
//                                          stat_length, hsp_list_out,
//                                          score_options->gapped_calculation,
//                                          isRPS, gap_align->sbp, 0, scale_factor);
//     }
//
//    /* Use score threshold rather than evalue if
//     * matrix_only_scoring is used.  -RMH-
//     */
//     if ( sbp->matrix_only_scoring )
//     {
//         status = Blast_HSPListReapByRawScore(hsp_list_out, hit_options);
//     }else {
//        /* Discard HSPs that don't pass the e-value test. */
//         status = s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
// ```
fn main() -> Result<()> {
    let cli: Cli = LOSAT::cli::try_parse_from(std::env::args_os())?;
    let Commands::Blastx(args) = cli.command else {
        bail!("BLASTX diagnostic only")
    };
    let options = args.resolve()?;
    let query = read_fasta(&args.query, false, args.lcase_masking)?;
    let subject = read_fasta(&args.subject, true, false)?;
    let db_length = subject.iter().map(|s| s.sequence.len() as i64).sum();
    for batch in search_preliminary(&query, &subject, &options)? {
        for chunk in &batch.chunks {
            for s in &chunk.subjects {
                let mut linked = evaluate_gapped(
                    &s.purged,
                    &chunk.chunk.prepared,
                    &s.parameters,
                    &options,
                    subject[s.oid].sequence.len() as i32,
                    db_length,
                    1.0,
                )?;
                dump("D_NUMERIC_PRE", s.call, s.oid, &linked);
                let deleted = reap_preliminary(&mut linked, &options);
                dump("D_RETAINED_PRE", s.call, s.oid, &linked);
                for h in deleted {
                    println!(
                        "D_DELETED_PRE\t{}\t{}\t{}\t{}\t{:016x}",
                        s.call,
                        s.oid,
                        h.source_index,
                        h.hsp.score,
                        h.evalue.to_bits()
                    );
                }
            }
        }
    }
    Ok(())
}
