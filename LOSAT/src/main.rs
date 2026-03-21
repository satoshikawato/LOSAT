use anyhow::Result;
use clap::{Parser, Subcommand};
use std::ffi::OsString;
use LOSAT::algorithm::{blastn, blastp, tblastx};

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

    /// Protein vs Protein
    Blastp(blastp::BlastpArgs),

    /// Translated DNA vs Translated DNA (Gapped)
    Tblastx(tblastx::TblastxArgs),
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
// ```c
// const string kArgQuery("query");
// const string kArgOutput("out");
// const string kArgSubject("subject");
// const string kTask("task");
// const string kArgNumThreads("num_threads");
// const string kArgMatrixName("matrix");
// const string kArgEvalue("evalue");
// const string kArgMaxTargetSequences("max_target_seqs");
// const string kArgGapOpen("gapopen");
// const string kArgGapExtend("gapextend");
// ```
fn normalize_ncbi_cli_args() -> Vec<OsString> {
    std::env::args_os()
        .map(|arg| {
            if let Some(text) = arg.to_str() {
                match text {
                    "-query" => OsString::from("-q"),
                    "-subject" => OsString::from("-s"),
                    "-out" => OsString::from("-o"),
                    "-outfmt" => OsString::from("--outfmt"),
                    "-evalue" => OsString::from("--evalue"),
                    "-task" => OsString::from("--task"),
                    "-num_threads" => OsString::from("-n"),
                    "-max_target_seqs" => OsString::from("--max-target-seqs"),
                    "-word_size" => OsString::from("--word-size"),
                    "-window_size" => OsString::from("--window-size"),
                    "-gapopen" => OsString::from("--gap-open"),
                    "-gapextend" => OsString::from("--gap-extend"),
                    "-comp_based_stats" => OsString::from("--comp-based-stats"),
                    "-matrix" => OsString::from("--matrix"),
                    "-seg" => OsString::from("--seg"),
                    _ => arg,
                }
            } else {
                arg
            }
        })
        .collect()
}

fn main() -> Result<()> {
    let startup_trace = std::env::var("LOSAT_STARTUP_TRACE").ok().as_deref() == Some("1");
    if startup_trace {
        eprintln!("[startup] enter main");
    }
    let cli = Cli::parse_from(normalize_ncbi_cli_args());
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
    }
    Ok(())
}
