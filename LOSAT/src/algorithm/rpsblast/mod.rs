//! RPSBLAST (Reverse Position-Specific BLAST): the user's protein queries against a database of
//! profiles (`makeprofiledb -dbtype rps`).
//!
//! LOSAT runs the search with and without composition-based statistics on databases of one or
//! several volumes, with blastp's driver (`blastp::args::BlastpApp::Rpsblast`: the query
//! batches and their masks, the database's statistics and reports) and the RPS engine
//! (`algorithm::rps_search`).

pub mod args;

pub use args::RpsblastArgs;

// NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:96-108
// ```c++
// int CRPSBlastApp::Run(void)
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
/// Runs the RPS-BLAST search with blastp's driver (`blastp::run`). LOSAT's output is NCBI's
/// `-num_threads 1` output at every thread count and `-mt_mode`.
pub fn run(
    args: RpsblastArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> anyhow::Result<()> {
    crate::algorithm::blastp::run(args.into_blastp_args()?, settings)
}
