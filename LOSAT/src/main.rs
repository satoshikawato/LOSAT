#![allow(warnings, clippy::all)]

use anyhow::Result;
use LOSAT::algorithm::{blastn, blastp, tblastn, tblastx};
use LOSAT::cli::{Cli, Commands};

fn main() -> Result<()> {
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
        Commands::Blastn(args) => {
            blastn::run(args)?;
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
