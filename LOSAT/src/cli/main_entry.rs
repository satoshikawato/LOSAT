//! The command line's entry: what the `losat` binary runs (`main_entry`). The binary names nothing
//! else of the library; the crate's root re-exports this function for it.

use std::io::Write;

use anyhow::Result;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
use crate::algorithm::{blastn, blastp, blastx, psiblast, tblastn, tblastx};
use crate::cli::{Cli, Commands};

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:407-410
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
/// where a panic ends the process it installs the panic hook of
/// `report::archive::strategy::install_panic_hook`, and on WASI it changes the working
/// directory), so a program that searches from its own code calls `api` instead.
#[doc(hidden)]
pub fn main_entry() -> Result<()> {
    // Where a panic (LOSAT's own failure) ends the process, it ends the strategy files as LOSAT's
    // refusals end them.
    crate::report::archive::strategy::install_panic_hook();
    let result = run_main_entry();
    // LOSAT's own errors (its refusals, which NCBI has no message for) stay on the terminal when
    // the standard error went to the `-logfile` target (tblastx, rpsblast and rpstblastn,
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
    // The process's run: NCBI's application starts with a diagnostic stream that has not failed,
    // the command line's (`ncbi_exception::main_diag_failbit`, `MainDiag`).
    // Diagnostic (LOSAT_TIMING): the clock of `since_start` (utils/xtiming.rs). No NCBI
    // counterpart; it reads the clock only with LOSAT_TIMING.
    crate::utils::xtiming::start();
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3719-3725
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
    // `losat blastdbcmd` reads its command line with NCBI's argument reader
    // (`ncbi_app::args`), not clap's: the words after `blastdbcmd` go to it as they are.
    // An executable called `blastn` is `losat blastn` (`cli::entry_argv`).
    let argv: Vec<std::ffi::OsString> = crate::cli::entry_argv(std::env::args_os());
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:84-86
    // Every NCBI BLAST+ program's application but the maskers' (dustmasker, segmasker,
    // windowmasker) has a usage
    // report, made with the application before `AppMain` runs; it reads `.ncbirc`
    // (`ncbi_environment::ApplicationStart`), and its messages go to stderr at once, before
    // `-version` and before blastdbcmd's `-logfile` is read. The start of the others (dustmasker,
    // segmasker, windowmasker, rpsbproc) writes nothing here.
    // The program's application start, made once here and given to what reads after it.
    let start = argv
        .get(1)
        .and_then(|word| word.to_str())
        .filter(|word| crate::cli::ENTRY_NAMES.contains(word))
        .map(|program| {
            crate::blastinput::ncbi_environment::ApplicationStart::of_process(program, &argv[2..])
        });
    // The start of the program `name`, which the entry made above when `argv` names it.
    let started = |name: &str| start.as_ref().filter(|start| start.program() == name);
    if let Some(start) = &start {
        match start.usage_report() {
            Ok(text) => {
                let mut stderr = std::io::stderr();
                let _ = stderr.write_all(text.as_bytes());
                let _ = stderr.flush();
            }
            Err(reason) => return Err(anyhow::Error::msg(reason.to_string())),
        }
    }
    if let Some(start) = started("blastdbcmd") {
        std::process::exit(crate::blastdbcmd::run(&argv[2..], start));
    }
    // NCBI reference (598d8ae6): c++/src/app/dustmasker/main.cpp:39-42
    // `losat dustmasker` reads its command line with NCBI's argument reader (`ncbi_app`), as
    // blastdbcmd does.
    if let Some(start) = started("dustmasker") {
        std::process::exit(crate::dustmasker::run(&argv[2..], start));
    }
    // NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:265-272
    // `losat segmasker` too.
    if let Some(start) = started("segmasker") {
        std::process::exit(crate::segmasker::run(&argv[2..], start));
    }
    // NCBI reference (598d8ae6): c++/src/app/winmasker/main.cpp:39-42
    // `losat windowmasker` too.
    if let Some(start) = started("windowmasker") {
        std::process::exit(crate::windowmasker::run(&argv[2..], start));
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:259-263
    // `losat convert2blastmask` too.
    if let Some(start) = started("convert2blastmask") {
        std::process::exit(crate::convert2blastmask::run(&argv[2..], start));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/main.cpp:51-54
    // `losat makembindex` too.
    if let Some(start) = started("makembindex") {
        std::process::exit(crate::makembindex::run(&argv[2..], start));
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1920-1923
    // `losat makeprofiledb` too.
    if let Some(start) = started("makeprofiledb") {
        std::process::exit(crate::makeprofiledb::run(&argv[2..], start));
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
    // The diagnostics go to the file of `-logfile`, which NCBI reads anywhere before `--`
    // (`ncbi_environment::log_file`), opened before `-version`: `blastn -version -logfile F`
    // makes F, with the usage report's messages.
    let mut target = match &start {
        Some(_) => {
            use crate::blastinput::ncbi_environment::{diag_target, log_file};
            diag_target(log_file(&argv[2..]))
        }
        None => crate::blastinput::ncbi_environment::DiagTarget::Stderr,
    };
    if let Some(start) = &start {
        start.write_log_setup(&mut target, &mut std::io::stderr());
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1001 (`-version` and the full
    // versions are read before the program's arguments and before the configuration load,
    // `cli::ncbi_version_request`).
    let argv_is_rpsbproc = argv.get(1).is_some_and(|word| word == "rpsbproc");
    if let Some(text) = crate::cli::ncbi_version_request(argv.iter().cloned()) {
        if argv_is_rpsbproc {
            crate::rpsbproc::write_stdout(&text);
        } else {
            write_dropping_errors(&mut std::io::stdout(), &text);
        }
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:413-424 (`x_TryInit`: the
    // configuration is loaded before the arguments are read; cited at
    // `ncbi_environment::ApplicationStart`). Its messages go where the diagnostics go, and a
    // `<program>.ini` that NCBI reports as malformed ends the program with exit status 2, `-h`
    // included.
    if let Some(start) = &start {
        match start.write_config_load(&mut target, &mut std::io::stderr()) {
            Ok(None) => {}
            Ok(Some(status)) => std::process::exit(status),
            Err(reason) => return Err(anyhow::Error::msg(reason)),
        }
    }
    // The top level (`losat`, `losat -h`, `losat --version`, a word that names no program) is
    // clap's.
    let Some(start) = start else {
        let top_level = crate::cli::top_level(&argv);
        match top_level.error {
            true => write_dropping_errors(&mut std::io::stderr(), &top_level.text),
            false => write_dropping_errors(&mut std::io::stdout(), &top_level.text),
        }
        std::process::exit(top_level.status);
    };
    let tool = crate::cli::reader::tool(start.program())
        .expect("a program of ENTRY_NAMES that is not on `ncbi_app` is on the reader");
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:685-697 (cited at `ncbi_app::read`)
    // NCBI's reader of the words that the application leaves to it (`-dryrun` taken), with what it
    // writes: the help on the standard output, an argument error's report (its usage on the
    // standard error, its posts where the diagnostics go), the XML of `-xmlhelp`
    // (`cli::xml_help_text`, without the arguments that LOSAT refuses).
    let (words, dryrun) = match crate::blastinput::ncbi_environment::preparse(&argv[2..]) {
        crate::blastinput::ncbi_environment::Preparsed::Words { words, dryrun, .. } => {
            (words, dryrun)
        }
        crate::blastinput::ncbi_environment::Preparsed::Version(_) => {
            unreachable!("the versions are printed above")
        }
    };
    let mut stderr = std::io::stderr();
    let mut diag = crate::ncbi_app::Diag::new(&target, &mut stderr, tool);
    let read = crate::ncbi_app::read(
        tool,
        &words,
        &mut std::io::stdout().lock(),
        &mut diag,
        Some(crate::cli::xml_help_text),
    );
    let args = match read {
        Ok(args) => args,
        Err(crate::ncbi_app::ReadEnd::Exit(status)) => std::process::exit(status),
        Err(crate::ncbi_app::ReadEnd::Refused(reason)) => std::process::exit(diag.refuse(&reason)),
    };
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:705-707, 292-296
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:424-431, 1035-1041
    // `-dryrun` (`ncbi_environment::preparse`): the arguments have passed their definitions (the
    // reader above), and NCBI's BLAST programs and makeblastdb keep the default `DryRun`: nothing
    // else runs, exit status 0 (NCBI BLAST+ 2.17.0: a missing `-query` or `-db` file, no
    // `-query`, `-out FILE` not made; line RG item E1b), before LOSAT's refusals (line RG
    // DECISIONS K4i-2b-2 (1)). Its `Info` post is written when the diagnostics' settings, applied
    // before it (`DIAG_POST_LEVEL`, the registry's `[DEBUG]`, `[Diag]`, `[LOG]`), let it through:
    // LOSAT refuses the settings it refuses in a run, which include those, with the same check
    // (line RG DECISIONS E1b-6). rpsbproc has its own `DryRun`, after its `Init`
    // (`rpsbproc::dry_run`, which makes the same check; item E1c), which reads its arguments.
    if dryrun && start.program() != "rpsbproc" {
        start.settings().map_err(anyhow::Error::msg)?;
        return Ok(());
    }
    // LOSAT's refusals of what it does not run (`reader::refusal_after_read`), then clap's
    // conversion of the arguments that NCBI's reader read (`reader::convert`): LOSAT's errors, on
    // the terminal's standard error (`Error: ` and the text, exit status 1; line RG DECISIONS 32).
    let command = match crate::cli::reader::refusal_after_read(tool, &args) {
        Some(refusal) => Err(refusal),
        None => crate::cli::reader::convert(tool, &args, false),
    };
    let command = match command {
        Ok(command) => command,
        Err(reason) if start.program() == "rpsbproc" => {
            exit_rpsbproc(crate::rpsbproc::Outcome::Refused(reason.into_bytes()));
            return Ok(());
        }
        Err(reason) => return Err(anyhow::Error::msg(reason)),
    };
    let cli = Cli { command };
    if dryrun {
        if let Commands::Rpsbproc(args) = cli.command {
            exit_rpsbproc(crate::rpsbproc::dry_run(
                args,
                start.settings(),
                crate::rpsbproc::Diagnostics::new(target),
            ));
        }
        return Ok(());
    }
    // The messages of tblastx's, rpsblast's and rpstblastn's runs go where the diagnostics go
    // (`ncbi_environment::run_messages_follow_log`): the standard error becomes the `-logfile`
    // target from here, where NCBI's `Run` makes its handler. LOSAT's refusals above stay on the
    // terminal's standard error (line RG DECISIONS 32). rpsbproc posts its
    // diagnostics to the target itself (`rpsbproc::Diagnostics`): its `cerr` lines stay on the
    // standard error.
    let mut rpsbproc_target = None;
    match start.program() {
        name if crate::blastinput::ncbi_environment::run_messages_follow_log(name) => {
            if let Err(reason) = target.into_standard_error(name) {
                return Err(anyhow::Error::msg(reason));
            }
        }
        "rpsbproc" => rpsbproc_target = Some(target),
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
        // Native serial search passes actual final results to the BLASTX formatter.
        Commands::Blastx(args) => {
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
            // NCBI's application layer reads its environment and registry files before blastx's
            // own code runs; LOSAT rejects the settings that change the output first, before any
            // option or file of the run, and `native::run` runs with them.
            let run = start
                .settings()
                .map_err(anyhow::Error::msg)
                .and_then(|settings| blastx::BlastxArgs::run(args, settings));
            if let Err(error) = run {
                blastx::native::exit_on_native_error(&error);
                return Err(error);
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-175,225-227
        Commands::Blastn(mut args) => {
            // The `@PG` line of SAM (outfmt 17) shows the command line.
            args.command_line = Some(crate::cli::command_line_bytes());
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1044
            // NCBI's application layer reads its environment and registry files before
            // blastn's own code runs; LOSAT rejects the settings that change the output.
            let settings = start.settings().map_err(anyhow::Error::msg)?;
            if let Err(error) = blastn::run(args, settings) {
                crate::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
        Commands::Blastp(args) => {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1044
            // NCBI's application layer reads its environment and registry files before
            // blastp's own code runs; LOSAT rejects the settings that change the output.
            let settings = start.settings().map_err(anyhow::Error::msg)?;
            if let Err(error) = blastp::run(args, settings) {
                crate::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227,252-255
        Commands::Tblastx(args) => {
            // NCBI's application layer reads its environment and registry files before
            // tblastx's own code runs (as for BLASTN); LOSAT rejects the settings that change
            // the output.
            let settings = start.settings().map_err(anyhow::Error::msg)?;
            if let Err(error) = tblastx::run(args, settings) {
                crate::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:288-301:
        // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1322-1333
        Commands::Makedb(args) => {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1041
            // makeblastdb is an NCBI application like the search programs: its application
            // layer reads the same environment and registry files (`makeblastdb.ini`, then
            // `.ncbirc`) before `Run`, and LOSAT rejects the settings that change the output.
            // The data loaders the check reports are not used for a FASTA input.
            let settings = start.settings().map_err(anyhow::Error::msg)?;
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
            // The variables that makeblastdb reads itself.
            crate::writedb::makeblastdb::check_environment(&args, &|name: &str| {
                std::env::var_os(name)
            })
            .map_err(anyhow::Error::msg)?;
            let status = crate::writedb::makeblastdb::run(&args, &settings);
            if status != 0 {
                std::process::exit(status);
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:583-585
        Commands::Psiblast(args) => {
            // psiblast is an NCBI application as blastp: its application layer reads its
            // environment and registry files (`psiblast.ini`, then `.ncbirc`) before its own
            // code runs; LOSAT rejects the settings that change the output.
            let settings = start.settings().map_err(anyhow::Error::msg)?;
            if let Err(error) = psiblast::run(args, settings) {
                crate::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:324-328
        Commands::Rpsblast(args) => {
            // rpsblast is an NCBI application as blastp: its application layer reads its
            // environment and registry files (`rpsblast.ini`, then `.ncbirc`) before its own
            // code runs; LOSAT rejects the settings that change the output.
            let settings = start.settings().map_err(anyhow::Error::msg)?;
            if let Err(error) = crate::algorithm::rpsblast::run(args, settings) {
                crate::cli::exit_on_native_error(&error);
                return Err(error);
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:319-322
        Commands::Rpstblastn(args) => {
            // rpstblastn is an NCBI application as blastx: its application layer reads its
            // environment and registry files (`rpstblastn.ini`, then `.ncbirc`) before its own
            // code runs; LOSAT rejects the settings that change the output, before its own
            // refusals. BLASTX's driver runs with them (`native::run`).
            let settings = start.settings().map_err(anyhow::Error::msg)?;
            if let Err(error) = crate::algorithm::rpstblastn::run(args, settings) {
                blastx::native::exit_on_native_error(&error);
                return Err(error);
            }
        }
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1870
        // rpsbproc ends with its own status (`exit_rpsbproc`).
        Commands::Rpsbproc(args) => {
            let target =
                rpsbproc_target.unwrap_or(crate::blastinput::ncbi_environment::DiagTarget::Stderr);
            exit_rpsbproc(crate::rpsbproc::run(
                args,
                start.settings(),
                crate::rpsbproc::Diagnostics::new(target),
            ));
        }
        Commands::Tblastn(args) => {
            // NCBI's application layer reads its environment and registry files before
            // tblastn's own code runs (as for BLASTN); LOSAT rejects the settings that change
            // the output.
            let settings = start.settings().map_err(anyhow::Error::msg)?;
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
    crate::algorithm::tblastx::lookup::compressed::x_clut_direct::x_clutpairs_print_summary();
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

/// Ends `losat rpsbproc` as its `outcome` says: its own status (`Run`'s or `DryRun`'s
/// `m_app_status`, or the application layer's for an exception; exit status 0 returns), or
/// LOSAT's refusal, an error (exit 1) written as `main` writes an error, through
/// `rpsbproc::write_stderr_bytes` (a closed standard error: exit 4; the message may quote data
/// bytes).
fn exit_rpsbproc(outcome: crate::rpsbproc::Outcome) {
    match outcome {
        crate::rpsbproc::Outcome::Exit(0) => {}
        crate::rpsbproc::Outcome::Exit(code) => std::process::exit(code),
        crate::rpsbproc::Outcome::Refused(message) => {
            let mut line = b"Error: ".to_vec();
            line.extend_from_slice(&message);
            line.push(b'\n');
            crate::rpsbproc::write_stderr_bytes(&line);
            std::process::exit(1);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:972-977 (cited at
// `ncbi_version::print`: the version is written to `cout`, whose failure NCBI never checks)
/// Writes `text` (a version, the top level's help or error) to `stream` as
/// NCBI writes its version, help and usage to `cout` and `cerr`: a write that fails is dropped
/// and the program ends with its exit status (NCBI BLAST+ 2.17.0 `blastn -version > /dev/full`,
/// `-version-full`, `-h`, `-help`: exit status 0 and no message, line RG item E3), where `print!`
/// panics (exit status 134). A pipe whose reader has closed ends NCBI's programs by SIGPIPE, which
/// LOSAT's runtime ignores (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exceptions 5 and 10).
fn write_dropping_errors(stream: &mut dyn std::io::Write, text: &str) {
    let _ = stream.write_all(text.as_bytes());
    let _ = stream.flush();
}
