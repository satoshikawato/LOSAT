//! DELTABLAST (Domain Enhanced Lookup Time Accelerated BLAST)
//!
//! deltablast is psiblast whose round 1 searches a PSSM made from the conserved domains that the
//! query hits in a domain database. LOSAT runs it as psiblast's rounds (`BlastpApp::Deltablast`,
//! `blastp::blast_engine::psi_round`): before round 1 of each query the domain search (the RPS
//! engine on `-rpsdb`, `blastp::blast_engine::delta_hooks`) finds the domains, `cdd_input` turns
//! its hits into the alignment of domains that the PSSM engine reads (`core::psi::cd_stages`), and
//! round 1 searches the PSSM it makes as a later round of psiblast does (`pssm_search`); the later
//! rounds are psiblast's.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:196-460

pub mod args;
pub(crate) mod cdd_input;

pub use args::DeltablastArgs;

/// deltablast's options besides psiblast's (`BlastpArgs::delta`).
#[derive(Debug, Clone, PartialEq)]
pub struct DeltaOptions {
    /// `-rpsdb`: the domain database.
    pub rpsdb: String,
    /// `-domain_inclusion_ethresh`: the E-value below which a domain hit joins the PSSM, and the
    /// E-value of the domain search.
    pub domain_inclusion_ethresh: f64,
    /// `-show_domain_hits` (refused by name, `ncbi_app::args::refusal`: item PS13b).
    pub show_domain_hits: bool,
}

// NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:196-460
/// Runs deltablast with blastp's engine (`blastp::run`).
pub fn run(
    args: DeltablastArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> anyhow::Result<()> {
    crate::algorithm::blastp::run(args.into_blastp_args()?, settings)
}
