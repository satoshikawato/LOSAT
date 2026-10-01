#![allow(warnings, clippy::all)]

use anyhow::Result;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
// ```c++
//     static const string kProgram("blastx");
//     arg.Reset(new CProgramDescriptionArgs(kProgram,
//                                   "Translated Query-Protein Subject BLAST"));
// ```
use LOSAT::algorithm::{blastn, blastp, blastx, tblastn, tblastx};
use LOSAT::cli::{Cli, Commands};

fn main() -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3719-3725
    // ```c++
    // bool CDir::SetCwd(const string& dir)
    // {
    //     if ( NcbiSys_chdir(_T_XCSTRING(dir)) != 0 ) {
    //         LOG_ERROR_ERRNO(51, "CDir::SetCwd(): Cannot change directory to: " + dir);
    //         return false;
    //     }
    //     return true;
    // ```
    // NCBI file readers resolve lexical filenames against the caller's cwd.
    // WASI starts at /; the existing command host supplies its cwd in PWD.
    // Set libc's cwd once before argument parsing, without rewriting filenames.
    #[cfg(target_os = "wasi")]
    if let Some(cwd) = std::env::var_os("PWD") {
        std::env::set_current_dir(&cwd)?;
    }
    let startup_trace = std::env::var("LOSAT_STARTUP_TRACE").ok().as_deref() == Some("1");
    if startup_trace {
        eprintln!("[startup] enter main");
    }
    // NCBI blastinput/cmdline_flags.cpp:46-94: kArgQuery("query"), kArgSubject("subject").
    let cli: Cli = match LOSAT::cli::try_parse_from(std::env::args_os()) {
        Ok(cli) => cli,
        Err(error) => {
            let message = LOSAT::cli::render_message(&error);
            if error.use_stderr() {
                eprint!("{message}");
            } else {
                print!("{message}");
            }
            std::process::exit(error.exit_code());
        }
    };
    if startup_trace {
        eprintln!("[startup] after clap parse");
    }

    match cli.command {
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-281
        // ```c++
        //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
        //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
        //                 results = lcl_blast.Run();
        // 	        BLAST_PROF_STOP( APP.LOOP.BLAST );
        //             }
        // ```
        // Native serial search passes actual final results to the BLASTX formatter.
        Commands::Blastx(args) => {
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
            // ```c++
            //     catch (const CArgException& e) {                                        \
            //         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
            //         exit_code = BLAST_INPUT_ERROR;                                      \
            //     }                                                                       \
            // ```
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
            // ```c++
            //         } else {                                                            \
            //             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
            //             exit_code = BLAST_ENGINE_ERROR;                                 \
            // ```
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
            // ```c++
            //     }                                                                       \
            //     catch (const std::ios::failure&) {                                      \
            //         LOG_POST(Error << "BLAST failed to write output");                  \
            //         exit_code = BLAST_OUTPUT_ERROR;                                     \
            //     }                                                                       \
            // ```
            if let Err(error) = blastx::BlastxArgs::run(args) {
                blastx::native::exit_on_native_error(&error);
                return Err(error);
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-175,225-227
        // ```c++
        //     catch (const blast::CInputException& e) {                               \
        //         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
        //         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
        //         exit_code = BLAST_INPUT_ERROR;                                      \
        // ...
        //             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
        //             exit_code = BLAST_ENGINE_ERROR;                                 \
        // ```
        Commands::Blastn(args) => {
            // NCBI reference: ncbi-blast/c++/src/corelib/ncbiapp.cpp:1031-1044
            // ```c
            //     // Setup some debugging features from environment variables.
            //     if ( !m_Environ->Get(DIAG_TRACE).empty() ) {
            //         SetDiagTrace(eDT_Enable, eDT_Enable);
            //     }
            //     string post_level = m_Environ->Get(DIAG_POST_LEVEL);
            // ```
            // NCBI's application layer reads its environment and registry files before
            // blastn's own code runs; LOSAT rejects the settings that change the output.
            LOSAT::blastinput::ncbi_environment::check_ncbi_application_settings("blastn")
                .map_err(anyhow::Error::msg)?;
            // NCBI's blastn keeps the C runtime's default action for SIGPIPE: a write to a
            // closed pipe ends it by the signal (oracle: exit 128 + 13 from the shell, no
            // message). The Rust runtime ignores SIGPIPE at startup; blastn restores it.
            #[cfg(unix)]
            restore_default_sigpipe();
            if let Err(error) = blastn::run(args) {
                LOSAT::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
        Commands::Blastp(args) => {
            blastp::run(args)?;
        }
        Commands::Tblastx(args) => {
            tblastx::run(args)?;
        }
        // NCBI c++/src/app/blast/tblastn_app.cpp:288-301:
        // results = lcl_blast.Run();
        // formatter.PrintOneResultSet(**result, query);
        Commands::Tblastn(args) => {
            tblastn::TblastnArgs::run(args)?;
        }
    }
    Ok(())
}

/// Restores the default action of SIGPIPE (the process ends when it writes to a closed
/// pipe), which the Rust runtime sets to "ignore" before `main`.
#[cfg(unix)]
fn restore_default_sigpipe() {
    extern "C" {
        fn signal(signum: i32, handler: usize) -> usize;
    }
    const SIGPIPE: i32 = 13;
    const SIG_DFL: usize = 0;
    // SAFETY: `signal` is the C library's function, which the Rust standard library
    // links on Unix; SIGPIPE is 13 and SIG_DFL is 0 on Linux and macOS. It changes only
    // the disposition of SIGPIPE, before any thread of the search starts.
    unsafe {
        signal(SIGPIPE, SIG_DFL);
    }
}
