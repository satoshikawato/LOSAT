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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-55
// ```c++
//     static const string kProgram("blastx");
//     arg.Reset(new CProgramDescriptionArgs(kProgram,
//                                   "Translated Query-Protein Subject BLAST"));
//     const bool kQueryIsProtein = false;
//     m_Args.push_back(arg);
//     m_ClientId = kProgram + " " + CBlastVersion().Print();
//
//     static const char kDefaultTask[] = "blastx";
//     SetTask(kDefaultTask);
// ```
pub mod blastx;
