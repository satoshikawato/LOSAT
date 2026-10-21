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
    // Diagnostic (LOSAT_TIMING): the clock of `since_start` (utils/xtiming.rs). No NCBI
    // counterpart; it reads the clock only with LOSAT_TIMING.
    LOSAT::utils::xtiming::start();
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
            let settings =
                LOSAT::blastinput::ncbi_environment::check_ncbi_application_settings("blastn")
                    .map_err(anyhow::Error::msg)?;
            if let Err(error) = blastn::run(args, settings) {
                LOSAT::cli::exit_on_native_error(&error);
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
                LOSAT::blastinput::ncbi_environment::check_ncbi_application_settings("blastp")
                    .map_err(anyhow::Error::msg)?;
            if let Err(error) = blastp::run(args, settings) {
                LOSAT::cli::exit_on_native_error(&error);
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
                LOSAT::blastinput::ncbi_environment::check_ncbi_application_settings("tblastx")
                    .map_err(anyhow::Error::msg)?;
            if let Err(error) = tblastx::run(args, settings) {
                LOSAT::cli::exit_on_native_error(&error);
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
            LOSAT::blastinput::ncbi_environment::check_ncbi_application_settings("makeblastdb")
                .map_err(anyhow::Error::msg)?;
            // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1134-1136
            // ```c++
            //     const string dont_scan_bioseq =
            //         GetEnvironment().Get("NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS");
            //     const bool scan_bioseq_4_cfastareader_usrobj = static_cast<bool>(dont_scan_bioseq.empty());
            // ```
            // The variables that makeblastdb reads itself.
            LOSAT::writedb::makeblastdb::check_environment(&|name: &str| std::env::var_os(name))
                .map_err(anyhow::Error::msg)?;
            let status = LOSAT::writedb::makeblastdb::run(&args);
            if status != 0 {
                std::process::exit(status);
            }
        }
        Commands::Tblastn(args) => {
            // NCBI's application layer reads its environment and registry files before
            // tblastn's own code runs (as for BLASTN); LOSAT rejects the settings that change
            // the output.
            let settings =
                LOSAT::blastinput::ncbi_environment::check_ncbi_application_settings("tblastn")
                    .map_err(anyhow::Error::msg)?;
            if let Err(error) = tblastn::TblastnArgs::run(args, settings) {
                LOSAT::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
    }
    LOSAT::utils::xstats::print();
    LOSAT::core::composition_adjustment::adjust_scores::x_adjmemo_print_stats();
    LOSAT::core::composition_adjustment::adjust_scores::x_newtonexact_print_stats();
    LOSAT::algorithm::tblastx::x_seed_bucket::print_shadow_stats();
    LOSAT::algorithm::tblastx::blast_engine::x_pair_par_print_summary();
    LOSAT::algorithm::tblastx::blast_engine::x_query_frames_print_summary();
    LOSAT::algorithm::tblastx::blast_engine::x_seed_release_print_summary();
    LOSAT::algorithm::tblastx::lookup::x_bplut_print_summary();
    LOSAT::algorithm::tblastx::lookup::x_tbxlut_print_summary();
    LOSAT::algorithm::tblastx::sum_stats_linking::x_hsp_inplace_print_summary();
    LOSAT::algorithm::tblastx::sum_stats_linking::x_link_sort_print_summary();
    LOSAT::algorithm::tblastx::sum_stats_linking::x_link_slices_print_summary();
    LOSAT::algorithm::tblastn::x_subject_translate_print_summary();
    LOSAT::algorithm::tblastn::x_subject_store_print_summary();
    LOSAT::algorithm::tblastn::x_scan_bf_print_summary();
    LOSAT::algorithm::tblastn::x_link_window_print_summary();
    LOSAT::algorithm::tblastn::x_tbn_hitbf_print_summary();
    LOSAT::algorithm::blastn::x_subject_store_print_summary();
    LOSAT::algorithm::blastn::x_query_dust_par_print_summary();
    LOSAT::algorithm::blastn::blast_engine::x_query_prep_par_print_summary();
    LOSAT::algorithm::blastn::blast_engine::x_qgroup_par_print_summary();
    LOSAT::algorithm::blastn::blast_engine::x_db_cut_print_summary();
    LOSAT::algorithm::blastp::x_bp_seed_print_summary();
    LOSAT::algorithm::tblastx::lookup::compressed::x_clut_direct::x_bpclut_print_summary();
    LOSAT::algorithm::tblastx::lookup::compressed::x_clut_direct::x_bplkpflat_print_summary();
    LOSAT::algorithm::blastp::x_bp_lkp_ahead_print_summary();
    LOSAT::algorithm::blastp::x_restricted_gap_print_summary();
    LOSAT::algorithm::blastn::alignment::greedy::x_greedy_fuse_print_summary();
    LOSAT::algorithm::blastn::lookup::x_disc_lookup_par::print_summary();
    LOSAT::algorithm::blastn::lookup::x_mb_lookup_print_summary();
    LOSAT::algorithm::blastp::x_bp_prelim_compact_print_summary();
    LOSAT::algorithm::blastp::x_bp_redo_par_print_summary();
    LOSAT::algorithm::blastp::x_bp_subj_run_print_summary();
    LOSAT::algorithm::blastp::x_bp_db_cut_print_summary();
    LOSAT::algorithm::blastp::x_bp_listed_lean_print_summary();
    LOSAT::algorithm::blastp::x_bp_redo_db_subj_print_summary();
    LOSAT::algorithm::blastp::x_bp_subj_meta_print_summary();
    LOSAT::algorithm::blastx::seed::x_seed_scan::print_shadow_stats();
    LOSAT::algorithm::blastp::x_ungapped_fast::print_shadow_stats();
    LOSAT::algorithm::tblastx::extension::x_tbxungap_print_shadow_stats();
    LOSAT::utils::x_logclone::print_shadow_stats();
    LOSAT::utils::x_seg_trim_memo::print_shadow_stats();
    LOSAT::seqdb::x_dbidxread_print_summary();
    // Diagnostic (LOSAT_TIMING): the common line, last (utils/xtiming.rs).
    LOSAT::utils::xtiming::print();
    Ok(())
}
