//! RPSTBLASTN (Translated Reverse Position-Specific BLAST): the six frames of the user's
//! nucleotide queries against a database of profiles (`makeprofiledb -dbtype rps`).
//!
//! LOSAT runs the search with BLASTX's driver and reports (`blastx::native`: the query batches,
//! their frames, strands, genetic codes and masks, the database's titles and every output format)
//! and the RPS engine (`algorithm::rps_search::translated`), which takes the place of BLASTX's own
//! search of the query batches (`blastx::args::RpsTblastnSearch`). Series PR, item PS12a
//! (DECISIONS 38 D1).

pub mod args;
pub(crate) mod search;

pub use args::RpstblastnArgs;

// NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:96-108
// ```c++
// int CRPSTBlastnApp::Run(void)
// {
// 	const CArgs& args = GetArgs();
// 	CMTArgs mt_args(args);
// 	if ((mt_args.GetMTMode() == CMTArgs::eSplitByQueries) &&
// 			(mt_args.GetNumThreads() > 1)){
// 		m_UsageReport.AddParam(CBlastUsageReport::eMTMode, args[kArgMTMode].AsInteger());
// 		return x_RunMTBySplitQuery();
// 	}
// 	else {
// 		return x_RunMTBySplitDB();
// 	}
// }
// ```
/// Runs the RPSTBLASTN search with BLASTX's driver (`blastx::native::run`) and the settings of
/// rpstblastn's application start. LOSAT's output is NCBI's `-num_threads 1` output at every
/// thread count and `-mt_mode`.
pub fn run(
    args: RpstblastnArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> anyhow::Result<()> {
    args.into_blastx_args()?.run(settings)
}
