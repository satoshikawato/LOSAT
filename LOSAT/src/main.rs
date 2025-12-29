use anyhow::Result;
use clap::{Parser, Subcommand};
use LOSAT::algorithm::{blastn, tblastx};

#[derive(Parser)]
#[command(name = "losat")]
#[command(version = "0.1.0")]
#[command(about = "A miniaturized reimplementation of BLAST algorithm", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {

    /// Nucleotide vs Nucleotide (Megablast/Blastn)
    Blastn(blastn::BlastnArgs),

    /// Translated DNA vs Translated DNA (Gapped)
    Tblastx(tblastx::TblastxArgs),

}

fn main() -> Result<()> {
    let startup_trace = std::env::var("LOSAT_STARTUP_TRACE").ok().as_deref() == Some("1");
    if startup_trace {
        eprintln!("[startup] enter main");
    }
    let cli = Cli::parse();
    if startup_trace {
        eprintln!("[startup] after clap parse");
    }

    match cli.command {
        Commands::Blastn(args) => {
            blastn::run(args)?;
        }
        Commands::Tblastx(args) => {
            tblastx::run(args)?;
        }
    }
    Ok(())
}
