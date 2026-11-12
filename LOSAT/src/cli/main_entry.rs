//! The command line's entry: what the `losat` binary runs (`main_entry`). The binary names nothing
//! else of the library; the crate's root re-exports this function for it.

use std::io::Write;

use anyhow::Result;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
// ```c++
//     static const string kProgram("blastx");
//     arg.Reset(new CProgramDescriptionArgs(kProgram,
//                                   "Translated Query-Protein Subject BLAST"));
// ```
use crate::algorithm::{blastn, blastp, blastx, psiblast, tblastn, tblastx};
use crate::cli::{Cli, Commands};

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:407-410
// ```c++
// int NcbiSys_main(int argc, ncbi::TXChar* argv[])
// {
//     return CBlastnApp().AppMain(argc, argv);
// }
// ```
// NCBI builds one executable per application, each a `main` that runs its application's
// `AppMain`. The `losat` binary is all of them (an executable named `blastn` is `losat blastn`,
// `cli::entry_argv`); the arms below cite each program's own entry.
/// The command `losat`: the body of the binary's `main` (`src/main.rs`), the only caller, which
/// returns what this returns. It reads the process's arguments and environment, runs the program
/// they name, and ends the way the command ends: `std::process::exit` with the program's status
/// where NCBI's application exits with one, or an `Err` that `main` returns (Rust then writes
/// `Error: ` and the error's `Debug` text to standard error and exits with status 1). It takes
/// nothing from its caller, so no caller can pass it a value it does not expect; it acts on the
/// whole process (it can end it, it sets `cli::post_errors_in_run` for the rest of the process,
/// and on WASI it changes the working directory), so a program that searches from its own code
/// calls `api` instead.
#[doc(hidden)]
pub fn main_entry() -> Result<()> {
    let result = run_main_entry();
    // LOSAT's own errors (its refusals, which NCBI has no message for) stay on the terminal when
    // the standard error went to the `-logfile` target (tblastx and rpsblast,
    // `ncbi_environment::run_messages_follow_log`): written as `main`'s `Result` writes them,
    // exit status 1. NCBI's errors are written by the run itself (`exit_on_native_error`), where
    // NCBI writes them. Line RG DECISIONS 32.
    if let Err(error) = &result {
        if let Some(mut terminal) = crate::blastinput::ncbi_environment::terminal_stderr() {
            let _ = writeln!(terminal, "Error: {error:?}");
            std::process::exit(1);
        }
    }
    result
}

fn run_main_entry() -> Result<()> {
    // Diagnostic (LOSAT_TIMING): the clock of `since_start` (utils/xtiming.rs). No NCBI
    // counterpart; it reads the clock only with LOSAT_TIMING.
    crate::utils::xtiming::start();
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
    // No NCBI counterpart: one stderr line per diagnostic variable that a build without the
    // diagnostics feature ignores (utils/xdiag.rs).
    for line in crate::utils::xdiag::ignored_notices() {
        eprintln!("{line}");
    }
    let startup_trace = crate::utils::xdiag::var("LOSAT_STARTUP_TRACE")
        .ok()
        .as_deref()
        == Some("1");
    if startup_trace {
        eprintln!("[startup] enter main");
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1295-1298
    // ```c++
    // int main(int argc, const char* argv[] /*, const char* envp[]*/)
    // {
    //     return CBlastDBCmdApp().AppMain(argc, argv);
    // }
    // ```
    // `losat blastdbcmd` reads its command line with NCBI's argument reader
    // (`ncbi_app::args`), not clap's: the words after `blastdbcmd` go to it as they are.
    // An executable called `blastn` is `losat blastn` (`cli::entry_argv`).
    let argv: Vec<std::ffi::OsString> = crate::cli::entry_argv(std::env::args_os());
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:84-86
    // ```c++
    // CBlastUsageReport::CBlastUsageReport()
    // {
    // 	x_CheckBlastUsageEnv();
    // ```
    // Every NCBI BLAST+ program's application but dustmasker's has a usage report, made with the
    // application before `AppMain` runs; it reads `.ncbirc` (`ncbi_environment::ApplicationStart`),
    // and its messages go to stderr at once, before `-version` and before blastdbcmd's `-logfile`
    // is read. The start of the others (dustmasker, rpsbproc) writes nothing here.
    let program = argv
        .get(1)
        .and_then(|word| word.to_str())
        .filter(|word| crate::cli::ENTRY_NAMES.contains(word));
    let start = program
        .map(|program| crate::blastinput::ncbi_environment::start_application(program, &argv[2..]));
    if let Some(start) = start {
        match start.usage_report() {
            Ok(text) => {
                let mut stderr = std::io::stderr();
                let _ = stderr.write_all(text.as_bytes());
                let _ = stderr.flush();
            }
            Err(reason) => return Err(anyhow::Error::msg(reason.to_string())),
        }
    }
    if argv.get(1).is_some_and(|word| word == "blastdbcmd") {
        std::process::exit(crate::blastdbcmd::run(&argv[2..]));
    }
    // NCBI reference (598d8ae6): c++/src/app/dustmasker/main.cpp:39-42
    // ```c++
    // int main( int argc, char * argv[] )
    // {
    //     return CDustMaskApplication().AppMain(argc, argv);
    // }
    // ```
    // `losat dustmasker` reads its command line with NCBI's argument reader (`ncbi_app`), as
    // blastdbcmd does.
    if argv.get(1).is_some_and(|word| word == "dustmasker") {
        std::process::exit(crate::dustmasker::run(&argv[2..]));
    }
    // No NCBI counterpart: makes those names (`cli::install_shims`). Read before clap, as the
    // top level has no options of its own.
    if argv
        .get(1)
        .is_some_and(|word| crate::cli::INSTALL_SHIMS_WORDS.iter().any(|w| word == w))
    {
        // An option word is not a directory (`./-h` names one).
        let option = argv
            .get(2)
            .is_some_and(|word| word.to_string_lossy().starts_with('-'));
        if argv.len() > 3 || option {
            eprintln!("usage: losat --install-shims [DIR]");
            std::process::exit(1);
        }
        match crate::cli::install_shims(argv.get(2).map(std::ffi::OsString::as_os_str)) {
            Ok(paths) => {
                for path in paths {
                    println!("{}", path.display());
                }
                return Ok(());
            }
            Err(message) => {
                eprintln!("losat --install-shims: {message}");
                std::process::exit(1);
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:868-871
    // ```c++
    //     // Setup logging as soon as possible.
    //     // Setup for diagnostics
    //     try {
    //         CDiagContext::SetupDiag(diag, 0, eDCM_NoChange, m_LogFile);
    // ```
    // The diagnostics go to the file of `-logfile`, which NCBI reads anywhere before `--`
    // (`ncbi_environment::log_file`), opened before `-version`: `blastn -version -logfile F`
    // makes F, with the usage report's messages. rpsbproc refuses it as it reads its arguments
    // (line RG item E1c), after the file is made, as NCBI makes it.
    let mut target = match program {
        Some(_) => {
            use crate::blastinput::ncbi_environment::{diag_target, log_file};
            diag_target(log_file(&argv[2..]))
        }
        None => crate::blastinput::ncbi_environment::DiagTarget::Stderr,
    };
    if let Some(start) = start {
        start.write_log_setup(&mut target);
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1001 (`-version` is read
    // before the program's arguments, `cli::ncbi_version_request`).
    let argv_is_rpsbproc = argv.get(1).is_some_and(|word| word == "rpsbproc");
    if let Some(text) = crate::cli::ncbi_version_request(argv.iter().cloned()) {
        if argv_is_rpsbproc {
            crate::rpsbproc::write_stdout(&text);
        } else {
            print!("{text}");
        }
        return Ok(());
    }
    // NCBI reads the other words of that scan there too (`ncbi_environment::preparse`): a full
    // version is printed, which LOSAT refuses (`cli::ncbi_toolkit_word_error`) before the
    // configuration load. The file of `-conffile` is that load's (`start_application`); rpsbproc
    // refuses it (line RG item E1c).
    if let Some(program) = program {
        use crate::blastinput::ncbi_environment::{preparse, Preparsed};
        let refused = match preparse(&argv[2..]) {
            Preparsed::Refused(word) => Some(word),
            Preparsed::Words {
                conffile: Some(_), ..
            } if program == "rpsbproc" => Some("-conffile"),
            _ => None,
        };
        if let Some(word) = refused {
            exit_on_parse_error(
                &crate::cli::ncbi_toolkit_word_error(program, word),
                argv_is_rpsbproc,
            );
        }
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:413-424 (`x_TryInit`: the
    // configuration is loaded before the arguments are read; quoted at
    // `ncbi_environment::ApplicationStart`). Its messages go where the diagnostics go, and a
    // `<program>.ini` that NCBI reports as malformed ends the program with exit status 2, `-h`
    // included.
    if let Some(start) = start {
        match start.write_config_load(&mut target, &mut std::io::stderr()) {
            Ok(None) => {}
            Ok(Some(status)) => std::process::exit(status),
            Err(reason) => return Err(anyhow::Error::msg(reason)),
        }
    }
    let program_name = program.map(str::to_string);
    let dryrun = program.is_some()
        && matches!(
            crate::blastinput::ncbi_environment::preparse(&argv[2..]),
            crate::blastinput::ncbi_environment::Preparsed::Words { dryrun: true, .. }
        );
    // NCBI blastinput/cmdline_flags.cpp:46-94: kArgQuery("query"), kArgSubject("subject").
    let cli: Cli = match crate::cli::try_parse_from(argv) {
        Ok(cli) => cli,
        Err(error) => exit_on_parse_error(&error, argv_is_rpsbproc),
    };
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:705-707, 292-296
    // ```c++
    //             try {
    //                 *exit_code = m_DryRun ? DryRun() : Run();
    //             }
    // ...
    // int CNcbiApplicationAPI::DryRun(void)
    // {
    //     ERR_POST_X(1, Info << "DryRun: default implementation does nothing");
    //     return 0;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:424-431, 1035-1041
    // ```c++
    //     CDiagContext::SetupDiag(diag, m_Config, eDCM_Flush, m_LogFile);
    //     CDiagContext::x_FinalizeSetupDiag();
    //
    //     // Setup the standard features from the config file.
    //     // Don't call till after LoadConfig()
    //     // NOTE: this will override environment variables,
    //     // except DIAG_POST_LEVEL which is Set*Fixed*.
    //     x_HonorStandardSettings();
    // ...
    //     string post_level = m_Environ->Get(DIAG_POST_LEVEL);
    //     if ( !post_level.empty() ) {
    //         EDiagSev sev;
    //         if (CNcbiDiag::StrToSeverityLevel(post_level.c_str(), sev)) {
    //             SetDiagFixedPostLevel(sev);
    //         }
    //     }
    // ```
    // `-dryrun` (`ncbi_environment::preparse`): the arguments have passed their definitions (the
    // parser above), and NCBI's BLAST programs, makeblastdb and blastdbcmd keep the default
    // `DryRun`: nothing else runs, exit status 0 (NCBI BLAST+ 2.17.0: a missing `-query` or `-db`
    // file, no `-query`, `-out FILE` not made; line RG item E1b). Its `Info` post is written when
    // the diagnostics' settings, applied before it (`DIAG_POST_LEVEL`, the registry's `[DEBUG]`,
    // `[Diag]`, `[LOG]`), let it through: LOSAT refuses the settings it refuses in a run, which
    // include those, with the same check (line RG DECISIONS E1b-6). rpsbproc refuses `-dryrun`
    // (item E1c).
    if dryrun {
        if let Some(name) = program_name.as_deref() {
            crate::blastinput::ncbi_environment::check_ncbi_application_settings(name)
                .map_err(anyhow::Error::msg)?;
        }
        return Ok(());
    }
    // The messages of tblastx's and rpsblast's runs go where the diagnostics go
    // (`ncbi_environment::run_messages_follow_log`): the standard error becomes the `-logfile`
    // target from here, where NCBI's `Run` makes its handler. LOSAT's parser errors above stay on
    // the standard error (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 1).
    match program_name.as_deref() {
        Some(name) if crate::blastinput::ncbi_environment::run_messages_follow_log(name) => {
            if let Err(reason) = target.into_standard_error(name) {
                return Err(anyhow::Error::msg(reason));
            }
        }
        _ => drop(target),
    }
    if startup_trace {
        eprintln!("[startup] after clap parse");
    }
    // The search programs write an error's text before its error archive, where NCBI's
    // `CATCH_ALL` posts it (`cli::post_native_error`).
    crate::cli::post_errors_in_run();

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
        Commands::Blastn(mut args) => {
            // The `@PG` line of SAM (outfmt 17) shows the command line.
            args.command_line = Some(crate::cli::command_line_bytes());
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
            let settings =
                crate::blastinput::ncbi_environment::check_ncbi_application_settings("blastn")
                    .map_err(anyhow::Error::msg)?;
            if let Err(error) = blastn::run(args, settings) {
                crate::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
        Commands::Blastp(args) => {
            // NCBI reference: ncbi-blast/c++/src/corelib/ncbiapp.cpp:1031-1044
            // ```c
            //     // Setup some debugging features from environment variables.
            //     if ( !m_Environ->Get(DIAG_TRACE).empty() ) {
            //         SetDiagTrace(eDT_Enable, eDT_Enable);
            //     }
            //     string post_level = m_Environ->Get(DIAG_POST_LEVEL);
            // ```
            // NCBI's application layer reads its environment and registry files before
            // blastp's own code runs; LOSAT rejects the settings that change the output.
            let settings =
                crate::blastinput::ncbi_environment::check_ncbi_application_settings("blastp")
                    .map_err(anyhow::Error::msg)?;
            if let Err(error) = blastp::run(args, settings) {
                crate::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
        // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:225-227,252-255
        // ```c
        //             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
        //             exit_code = BLAST_ENGINE_ERROR;                                 \
        //     ...
        //     catch (const std::ios::failure&) {                                      \
        //         LOG_POST(Error << "BLAST failed to write output");                  \
        //         exit_code = BLAST_OUTPUT_ERROR;                                     \
        // ```
        Commands::Tblastx(args) => {
            // NCBI's application layer reads its environment and registry files before
            // tblastx's own code runs (as for BLASTN); LOSAT rejects the settings that change
            // the output.
            let settings =
                crate::blastinput::ncbi_environment::check_ncbi_application_settings("tblastx")
                    .map_err(anyhow::Error::msg)?;
            if let Err(error) = tblastx::run(args, settings) {
                crate::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
        // NCBI c++/src/app/blast/tblastn_app.cpp:288-301:
        // results = lcl_blast.Run();
        // formatter.PrintOneResultSet(**result, query);
        // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1322-1333
        // ```c++
        // int CMakeBlastDBApp::Run(void)
        // {
        //     SetDiagPostLevel(eDiag_Warning);
        //     SetDiagPostPrefix("makeblastdb");
        //
        //     int status = 0;
        //     try { x_BuildDatabase(); }
        //     CATCH_ALL(status)
        //     x_AddCmdOptions();
        //     m_UsageReport.AddParam(CBlastUsageReport::eExitStatus, status);
        //     return status;
        // }
        // ```
        Commands::Makedb(args) => {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1041
            // ```c++
            //     // Setup some debugging features from environment variables.
            //     if ( !m_Environ->Get(DIAG_TRACE).empty() ) {
            //         SetDiagTrace(eDT_Enable, eDT_Enable);
            //     }
            //     string post_level = m_Environ->Get(DIAG_POST_LEVEL);
            //     if ( !post_level.empty() ) {
            //         EDiagSev sev;
            //         if (CNcbiDiag::StrToSeverityLevel(post_level.c_str(), sev)) {
            //             SetDiagFixedPostLevel(sev);
            //         }
            //     }
            // ```
            // makeblastdb is an NCBI application like the search programs: its application
            // layer reads the same environment and registry files (`makeblastdb.ini`, then
            // `.ncbirc`) before `Run`, and LOSAT rejects the settings that change the output.
            // The data loaders the check reports are not used for a FASTA input.
            let settings =
                crate::blastinput::ncbi_environment::check_ncbi_application_settings("makeblastdb")
                    .map_err(anyhow::Error::msg)?;
            crate::writedb::makeblastdb::check_registry_settings(
                &args,
                settings
                    .max_seqid_length
                    .as_ref()
                    .map(Option::as_deref)
                    .map_err(String::as_str),
            )
            .map_err(anyhow::Error::msg)?;
            // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1134-1136
            // ```c++
            //     const string dont_scan_bioseq =
            //         GetEnvironment().Get("NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS");
            //     const bool scan_bioseq_4_cfastareader_usrobj = static_cast<bool>(dont_scan_bioseq.empty());
            // ```
            // The variables that makeblastdb reads itself.
            crate::writedb::makeblastdb::check_environment(&|name: &str| std::env::var_os(name))
                .map_err(anyhow::Error::msg)?;
            let status = crate::writedb::makeblastdb::run(&args);
            if status != 0 {
                std::process::exit(status);
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:583-585
        // ```c++
        // int main(int argc, const char* argv[] /*, const char* envp[]*/)
        // {
        //     return CPsiBlastApp().AppMain(argc, argv);
        // ```
        Commands::Psiblast(args) => {
            // psiblast is an NCBI application as blastp: its application layer reads its
            // environment and registry files (`psiblast.ini`, then `.ncbirc`) before its own
            // code runs; LOSAT rejects the settings that change the output.
            let settings =
                crate::blastinput::ncbi_environment::check_ncbi_application_settings("psiblast")
                    .map_err(anyhow::Error::msg)?;
            if let Err(error) = psiblast::run(args, settings) {
                crate::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:324-328
        // ```c++
        // int main(int argc, const char* argv[] /*, const char* envp[]*/)
        // {
        // 	int status = CRPSBlastApp().AppMain(argc, argv);
        // ```
        Commands::Rpsblast(args) => {
            // rpsblast is an NCBI application as blastp: its application layer reads its
            // environment and registry files (`rpsblast.ini`, then `.ncbirc`) before its own
            // code runs; LOSAT rejects the settings that change the output.
            let settings =
                crate::blastinput::ncbi_environment::check_ncbi_application_settings("rpsblast")
                    .map_err(anyhow::Error::msg)?;
            if let Err(error) = crate::algorithm::rpsblast::run(args, settings) {
                crate::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1870
        // ```c++
        // 	return CRpsbProcApp().AppMain(argc, argv);
        // ```
        // rpsbproc ends with its own status (`Run`'s `m_app_status`, or the application layer's
        // for an exception); LOSAT's refusals are errors (exit 1), written as `main` writes an
        // error, through `rpsbproc::write_stderr_bytes` (a closed standard error: exit 4; the
        // message may quote data bytes).
        Commands::Rpsbproc(args) => match crate::rpsbproc::run(args) {
            crate::rpsbproc::Outcome::Exit(0) => {}
            crate::rpsbproc::Outcome::Exit(code) => std::process::exit(code),
            crate::rpsbproc::Outcome::Refused(message) => {
                let mut line = b"Error: ".to_vec();
                line.extend_from_slice(&message);
                line.push(b'\n');
                crate::rpsbproc::write_stderr_bytes(&line);
                std::process::exit(1);
            }
        },
        Commands::Tblastn(args) => {
            // NCBI's application layer reads its environment and registry files before
            // tblastn's own code runs (as for BLASTN); LOSAT rejects the settings that change
            // the output.
            let settings =
                crate::blastinput::ncbi_environment::check_ncbi_application_settings("tblastn")
                    .map_err(anyhow::Error::msg)?;
            if let Err(error) = tblastn::TblastnArgs::run(args, settings) {
                crate::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
    }
    // Diagnostic (LOSAT_X_REACH): the sides the conditional switches took (utils/xreach.rs). No
    // NCBI counterpart; it writes only with LOSAT_X_REACH.
    crate::utils::xreach::write_report();
    crate::utils::xstats::print();
    crate::core::composition_adjustment::adjust_scores::x_adjmemo_print_stats();
    crate::core::composition_adjustment::adjust_scores::x_newtonexact_print_stats();
    crate::algorithm::tblastx::x_seed_bucket::print_summary();
    crate::algorithm::tblastx::blast_engine::x_pair_par_print_summary();
    crate::algorithm::tblastx::blast_engine::x_seed_release_print_summary();
    crate::algorithm::tblastx::sum_stats_linking::x_hsp_inplace_print_summary();
    crate::algorithm::tblastx::sum_stats_linking::x_link_sort_print_summary();
    crate::algorithm::tblastx::sum_stats_linking::x_link_slices_print_summary();
    crate::algorithm::tblastn::x_subject_translate_print_summary();
    crate::algorithm::tblastn::x_subject_store_print_summary();
    crate::algorithm::tblastn::x_scan_bf_print_summary();
    crate::algorithm::tblastn::x_scan_byte_print_summary();
    crate::algorithm::tblastn::x_lut_lean_print_summary();
    crate::algorithm::tblastn::x_link_window_print_summary();
    crate::algorithm::tblastn::x_link_ctx_print_summary();
    crate::algorithm::blastn::x_subject_store_print_summary();
    crate::algorithm::blastn::x_query_dust_par_print_summary();
    crate::algorithm::blastn::x_greedy_row_print_summary();
    crate::algorithm::blastn::blast_engine::x_query_prep_par_print_summary();
    crate::algorithm::blastn::blast_engine::x_qgroup_par_print_summary();
    crate::algorithm::blastn::blast_engine::x_tb_range_print_summary();
    crate::algorithm::blastn::blast_engine::x_lut_drop_print_summary();
    crate::algorithm::blastn::blast_engine::x_query_drop_print_summary();
    crate::algorithm::blastp::x_bp_hit_bf_print_summary();
    crate::algorithm::blastp::x_bp_flat_print_summary();
    crate::algorithm::tblastx::lookup::compressed::x_clut_direct::x_bplkpflat_print_summary();
    crate::algorithm::blastp::x_bp_lkp_ahead_print_summary();
    crate::algorithm::blastp::x_restricted_gap_print_summary();
    crate::algorithm::blastn::lookup::x_disc_lookup_par::print_summary();
    crate::algorithm::blastn::lookup::x_dc_csr::print_summary();
    crate::algorithm::blastn::lookup::x_mb_lookup_print_summary();
    crate::algorithm::blastx::seed::x_seed_scan::print_shadow_stats();
    crate::utils::x_logclone::print_shadow_stats();
    // Diagnostic (LOSAT_TIMING): the common line, last (utils/xtiming.rs).
    crate::utils::xtiming::print();
    Ok(())
}

/// No NCBI counterpart: an error of LOSAT's argument parser (`cli::try_parse_from`) or one of its
/// refusals (`cli::ncbi_toolkit_word_error`), written as clap writes it, and its exit status.
fn exit_on_parse_error(error: &clap::Error, argv_is_rpsbproc: bool) -> ! {
    // rpsbproc's arguments are not the single-dash words of the BLAST programs.
    let message = if argv_is_rpsbproc {
        crate::cli::render_rpsbproc_message(error)
    } else {
        crate::cli::render_message(error)
    };
    // rpsbproc's messages and help are written without a panic on a closed stream
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 9).
    if argv_is_rpsbproc {
        if error.use_stderr() {
            crate::rpsbproc::write_stderr(&message);
        } else {
            crate::rpsbproc::write_stdout(&message);
        }
    } else if error.use_stderr() {
        eprint!("{message}");
    } else {
        print!("{message}");
    }
    std::process::exit(error.exit_code());
}
