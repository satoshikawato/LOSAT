// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blastp_args.cpp:71-115
// ```c
// arg.Reset(new CGenericSearchArgs(kQueryIsProtein, false, false,
//                                  false, false, true));
// ...
// arg.Reset(new CWindowSizeArg);
// ...
// arg.Reset(new CCompositionBasedStatsArgs);
// ```
pub mod blastn;
pub mod blastp;
pub mod common;
pub mod tblastx;
