//! PSIBLAST (Position-Specific Iterated BLAST)
//!
//! LOSAT runs psiblast's first round as the search of each query with blastp's engine, as NCBI's
//! `CPsiBlast` runs it before any PSSM exists (`ePSIBlast`, one query per batch), with
//! psiblast's report (`blastp::args::BlastpApp::Psiblast`); with `-num_iterations` other than 1
//! each later round searches the PSSM that the engine made from the round before (`rounds`,
//! `blastp::blast_engine::psi_round`); `-in_pssm` searches a PSSM query (`in_pssm`).

pub mod args;
pub(crate) mod in_pssm;
pub(crate) mod in_pssm_files;
pub mod pssm;
pub mod rounds;

pub use args::PsiblastArgs;

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:233-241
// ```c++
//     CRef<IQueryFactory> query_factory(s_MakeQueryFactory(pssm, query));
//             CRef<CPsiBlast> psiblast;
//             if (run_token & 2)
//             {
//                 _ASSERT(psi_opts.NotEmpty());
//                 if (pssm.NotEmpty()) {
//                     psiblast.Reset(new CPsiBlast(pssm, db_adapter, psi_opts));
//                 } else {
//                     psiblast.Reset(new CPsiBlast(query_factory, db_adapter, psi_opts));
// ```
/// Runs psiblast's first round with blastp's engine (`blastp::run`).
pub fn run(
    args: PsiblastArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> anyhow::Result<()> {
    crate::algorithm::blastp::run(args.into_blastp_args()?, settings)
}
