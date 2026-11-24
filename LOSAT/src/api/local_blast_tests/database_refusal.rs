//! Every program's `run_local` refuses `-db` with one text (`RUN_LOCAL_DB_REFUSAL`) before it
//! writes anything: `run_local` searches the caller's subject records, and a `-db` search runs
//! through the program's `run_local_database`.
//!
//! No NCBI counterpart: the caller contract of LOSAT's entries.

use crate::api::local_blast::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_tblastn, run_local_tblastx,
    OutputSink, ReportOutputs,
};
use crate::api::QueryEnd;
use crate::cli::{try_parse_from, Cli, Commands};

#[test]
fn every_program_refuses_db_in_run_local_with_one_text() {
    for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
        for outfmt in ["6", "11"] {
            let label = format!("{program} outfmt {outfmt}");
            let cli: Cli = try_parse_from(["losat", program, "-query", "q", "-db", "nt"])
                .unwrap_or_else(|error| panic!("{label}: {error}"));
            let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
            let result = {
                let mut outputs =
                    ReportOutputs::single(outfmt, OutputSink::Writer(&mut sink), &mut diagnostics);
                let end = QueryEnd::Input;
                match cli.command {
                    Commands::Blastn(args) => run_local_blastn(args, &[], &end, &[], &mut outputs),
                    Commands::Blastp(args) => {
                        run_local_blastp(args, &[], &end, &[], "", "", &mut outputs)
                    }
                    Commands::Blastx(args) => run_local_blastx(args, &[], &end, &[], &mut outputs),
                    Commands::Tblastn(args) => {
                        run_local_tblastn(args, &[], &end, &[], &mut outputs)
                    }
                    Commands::Tblastx(args) => {
                        run_local_tblastx(args, &[], &end, &[], &mut outputs)
                    }
                    _ => unreachable!("a search program"),
                }
            };
            assert_eq!(
                format!("{:#}", result.expect_err(&label)),
                "-db (a database search) is run by run_local_database, not run_local",
                "{label}"
            );
            assert!(sink.is_empty() && diagnostics.is_empty(), "{label}");
        }
    }
}
