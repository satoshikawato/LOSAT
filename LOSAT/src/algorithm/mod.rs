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
// NCBI reference: c++/src/algo/blast/blastinput/tblastn_args.cpp:45-62
// static const string kProgram("tblastn");
// static const char kDefaultTask[] = "tblastn";
// SetTask(kDefaultTask);
pub mod tblastn;
pub mod tblastx;
