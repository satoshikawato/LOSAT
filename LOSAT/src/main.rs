use anyhow::Result;
use clap::{Parser, Subcommand};
use rust_blast::algorithm::{blastn, tblastx};

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
    let cli = Cli::parse();

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
