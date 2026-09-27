//! Comparison-only native caller; final acceptance uses the public LOSAT binary.
use anyhow::{bail, Result};
use LOSAT::cli::{Cli, Commands};
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:271-297
// ```c++
//                     InitializeRemoteBlast(queries, db_args, m_OptsHndl,
//                           m_CmdLineArgs->ProduceDebugRemoteOutput(),
//                           m_CmdLineArgs->GetClientId());
//                 results = rmt_blast->GetResultSet();
//             } else {
// 	        BLAST_PROF_START( APP.LOOP.BLAST );
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
// 	        BLAST_PROF_STOP( APP.LOOP.BLAST );
//             }
// 	    BLAST_PROF_START( APP.LOOP.FMT );
//             if (fmt_args->ArchiveFormatRequested(args)) {
//                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
//                 m_Bah.ResetMessages();
//             } else {
//                 BlastFormatter_PreFetchSequenceData(*results, scope,
//                 		                            fmt_args->GetFormattedOutputChoice());
//             	ITERATE(CSearchResultSet, result, *results) {
//                	    formatter.PrintOneResultSet(**result, query_batch);
//             	}
//             }
// 	    BLAST_PROF_STOP( APP.LOOP.FMT );
// 	    batch_num++;
//         }
//         BLAST_PROF_START( APP.POST );
//         formatter.PrintEpilog(opt);
// ```
fn main() -> Result<()> {
    let cli: Cli = LOSAT::cli::try_parse_from(std::env::args_os())?;
    let Commands::Blastx(args) = cli.command else {
        bail!("requires blastx");
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-184
    // ```c++
    //     catch (const CArgException& e) {                                        \
    //         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
    //         exit_code = BLAST_INPUT_ERROR;                                      \
    //     }                                                                       \
    //     catch (const CObjReaderParseException& e) {                             \
    //         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
    //         exit_code = BLAST_INPUT_ERROR;                                      \
    //     }                                                                       \
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
    // ```c++
    //         } else {                                                            \
    //             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
    //             exit_code = BLAST_ENGINE_ERROR;                                 \
    // ```
    if let Err(error) = LOSAT::algorithm::blastx::native::run(&args) {
        LOSAT::algorithm::blastx::native::exit_on_native_error(&error);
        return Err(error);
    }
    Ok(())
}
