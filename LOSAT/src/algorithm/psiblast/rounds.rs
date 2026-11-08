//! psiblast's rounds of a FASTA query (`-num_iterations` other than 1; series PS, item PS5a):
//! the plan of a run (`PsiRoundPlan`), NCBI's iteration state (`IterationState`, a port of
//! `CPsiBlastIterationState`), the inclusion set of a round (`InclusionSet`), what one round's
//! search hands the round driver (`PsiRoundIo`, `RoundOutcome`) and the PSSM files the driver
//! writes between rounds (`RoundFiles`). The driver itself runs blastp's engine once per round per
//! query (`blastp::blast_engine::psi_round`).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-321
//! ```c++
//!             while (itr) {
//!
//!                 SavePssmToFile(pssm, &itr);
//!
//!                 CRef<CSearchResultSet> results;
//!                 ...
//!                    results = psiblast->Run();
//!                    ...
//! 	                      formatter.PrintOneResultSet(**result, query,
//! 	                                                itr.GetIterationNumber(),
//! 	                                                itr.GetPreviouslyFoundSeqIds());
//!                 ...
//!                 CSearchResults& results_1st_query = (*results)[0];
//!                 if ( !results_1st_query.HasAlignments() ) {
//!                     break;
//!                 }
//!
//!                 if (run_token & 2)
//!                 {
//!                         CConstRef<CSeq_align_set> aln(results_1st_query.GetSeqAlign());
//!                         CPsiBlastIterationState::TSeqIds ids;
//!                         CPsiBlastIterationState::GetSeqIds(aln, psi_opts, ids);
//!
//!                         itr.Advance(ids);
//!
//!                         if (itr || m_CmdLineArgs->GetSaveLastPssm()) {
//!                         CConstRef<CBioseq> seq =
//!                          s_GetQueryBioseq(query, scope, pssm);
//!                         pssm =
//!                          ComputePssmForNextIteration(*seq, aln, psi_opts, scope,
//!                                      results_1st_query.GetAncillaryData());
//!                         psiblast->SetPssm(pssm);
//!                         }
//!                 }
//!                 else
//!                   break;
//!             }
//!             if (itr.HasConverged())
//!               converged = true;
//!
//!        if (m_CmdLineArgs->GetSaveLastPssm()) {
//!            SavePssmToFile(pssm, &itr);
//!        }
//! ```

use std::collections::BTreeSet;
use std::io::Write as _;

use anyhow::Result;

use super::pssm::{
    is_regular_file, OutputFileReset, PsiblastPssmOptions, PssmFileCheck, PssmRequest, QueryPssm,
};

/// The rounds of a psiblast run of FASTA queries: `-num_iterations` (0: until convergence), the
/// PSSM engine's options and files (`-pseudocount`, `-out_pssm`, `-out_ascii_pssm`,
/// `-save_each_pssm`) and `-save_pssm_after_last_round`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PsiRoundPlan {
    pub(crate) num_iterations: u32,
    pub(crate) request: PssmRequest,
    pub(crate) save_last: bool,
}

impl PsiRoundPlan {
    /// The plan of psiblast's options with `-num_iterations` other than 1 (`None` with 1).
    pub(crate) fn of(options: &PsiblastPssmOptions) -> Option<Self> {
        let num_iterations = options.rounds?;
        Some(Self {
            num_iterations,
            request: PssmRequest {
                asn_path: options.out_pssm.clone(),
                ascii_path: options.out_ascii_pssm.clone(),
                save_each: options.save_each_pssm,
                pseudocount: options.pseudocount,
            },
            save_last: options.save_pssm_after_last_round,
        })
    }

    /// The PSSM files of the run, for the checks before the search (`check_pssm_titles`,
    /// `check_pssm_files`); `None` without `-out_pssm` and `-out_ascii_pssm`.
    pub(crate) fn file_request(&self) -> Option<&PssmRequest> {
        (self.request.asn_path.is_some() || self.request.ascii_path.is_some())
            .then_some(&self.request)
    }
}

/// The key of a subject in the inclusion sets (`CSeq_id_Handle` of its Seq-align id, DESIGN_ps5
/// §2): the index of a `-subject` record (`lcl|Subject_N`), the ordinal id of a database sequence
/// without parsed Seq-ids (`gnl|BL_ORD_ID|oid`), or the FASTA string of a parsed Seq-align id
/// (copies of one record that the reports join share it).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum SetKey {
    Index(u32),
    Id(Vec<u8>),
}

/// A round's inclusion set (`CPsiBlastIterationState::TSeqIds`).
pub(crate) type InclusionSet = BTreeSet<SetKey>;

/// `CPsiBlastIterationState`.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_iteration.cpp:49-130
// ```c++
// CPsiBlastIterationState::CPsiBlastIterationState(unsigned int num_iterations)
//     : m_TotalNumIterationsToDo(num_iterations), m_IterationsDone(0)
// {}
// ...
// bool
// CPsiBlastIterationState::HasMoreIterations() const
// {
//     if (m_TotalNumIterationsToDo == 0 ||
//         m_IterationsDone < m_TotalNumIterationsToDo) {
//         return true;
//     }
//     return false;
// }
//
// bool
// CPsiBlastIterationState::HasConverged()
// {
//     // For an object that hasn't been 'advanced' or one that only has performed
//     // one iteration (i.e.: called Advance() once), it doesn't make sense to
//     // have converged
//     if (m_IterationsDone <= 1) {
//         return false;
//     }
//     // if the most current list of ids is empty, that means that no new
//     // sequences were found and therefore we have converged.
//     if ( !m_PreviousData.empty() && m_CurrentData.empty() ) {
//         return true;
//     }
//
//     // Element by element comparison: has every element in the current Seq-id
//     // list been found in the previous iterations?
//     bool retval = true;
//     ITERATE(TSeqIds, current_id, m_CurrentData) {
//         if (m_PreviousData.find(*current_id) == m_PreviousData.end()) {
//             retval = false;
//             break;
//         }
//     }
//     return retval;
// }
//
// CPsiBlastIterationState::operator bool()
// {
//     return (HasMoreIterations() && !HasConverged());
// }
// ...
// CPsiBlastIterationState::TSeqIds
// CPsiBlastIterationState::GetPreviouslyFoundSeqIds() const
// {
//     return m_CurrentData;
// }
//
// void
// CPsiBlastIterationState::Advance(const TSeqIds& list)
// {
//     x_ThrowExceptionOnLogicError();
//
//     m_PreviousData = m_CurrentData;
//     m_CurrentData = list;
//     m_IterationsDone++;
// }
//
// unsigned int
// CPsiBlastIterationState::GetIterationNumber() const
// {
//     return m_IterationsDone+1;
// }
// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct IterationState {
    total: u32,
    done: u32,
    previous: InclusionSet,
    current: InclusionSet,
}

impl IterationState {
    pub(crate) fn new(num_iterations: u32) -> Self {
        Self {
            total: num_iterations,
            ..Self::default()
        }
    }

    pub(crate) fn has_more_iterations(&self) -> bool {
        self.total == 0 || self.done < self.total
    }

    pub(crate) fn has_converged(&self) -> bool {
        if self.done <= 1 {
            return false;
        }
        if !self.previous.is_empty() && self.current.is_empty() {
            return true;
        }
        self.current.is_subset(&self.previous)
    }

    /// `operator bool`: another round follows.
    pub(crate) fn go(&self) -> bool {
        self.has_more_iterations() && !self.has_converged()
    }

    /// `Advance` (`x_ThrowExceptionOnLogicError` is never reached: the driver advances only
    /// while `go`).
    pub(crate) fn advance(&mut self, set: InclusionSet) {
        debug_assert!(self.go());
        self.previous = std::mem::replace(&mut self.current, set);
        self.done += 1;
    }

    /// `GetIterationNumber`: the number of the round to run.
    pub(crate) fn number(&self) -> u32 {
        self.done + 1
    }

    /// `GetPreviouslyFoundSeqIds`: the inclusion set of the round before the current one.
    pub(crate) fn previously_found(&self) -> &InclusionSet {
        &self.current
    }
}

/// What the round driver hands one round's search (`BlastpReports::psi_round`), and what the
/// search hands back (`outcome`, filled by `blast_engine::psi_round::round_outcome` from the
/// round's Seq-aligns before the formatter's cut).
pub(crate) struct PsiRoundIo {
    /// The iteration state before this round's `Advance` (its number and previous set).
    pub(crate) state: IterationState,
    pub(crate) save_last: bool,
    /// The PSSM engine's options and the files' requests (`PsiRoundPlan::request`).
    pub(crate) request: PssmRequest,
    pub(crate) outcome: std::sync::Mutex<Option<RoundOutcome>>,
}

impl PsiRoundIo {
    /// Whether the round makes the PSSM of the next round, or the last one that is saved
    /// (`RoundOutcome::next`): it reads the letters of every subject of the round's Seq-aligns
    /// before the formatter's cut, by its Seq-align id.
    pub(crate) fn makes_next(&self) -> bool {
        self.outcome
            .lock()
            .expect("the round's outcome")
            .as_ref()
            .is_some_and(|outcome| outcome.next.is_some())
    }
}

/// The outcome of a round: whether it has alignments (`HasAlignments`), its inclusion set
/// (`GetSeqIds`), and the PSSM made from its Seq-aligns where NCBI makes one (after `Advance`,
/// `itr || save_last`), with NCBI's error of the PSSM engine, which comes after the round's
/// report.
pub(crate) struct RoundOutcome {
    pub(crate) has_alignments: bool,
    pub(crate) set: InclusionSet,
    pub(crate) next: Option<Result<QueryPssm>>,
}

/// The PSSM files of a run (`SavePssmToFile`), written between the rounds: each option's
/// `CAutoOutputFileReset` (one run-level `F.N` counter each), the ASN.1 file then the ASCII one.
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:385-407 (quoted at
// `pssm::PssmFiles`)
pub(crate) struct RoundFiles {
    asn: Option<OutputFileReset>,
    ascii: Option<OutputFileReset>,
    check: PssmFileCheck,
}

impl RoundFiles {
    pub(crate) fn new(request: &PssmRequest, out: Option<&std::path::Path>) -> Self {
        let stream =
            |path: Option<&str>| path.map(|path| OutputFileReset::new(path, request.save_each));
        Self {
            asn: stream(request.asn_path.as_deref()),
            ascii: stream(request.ascii_path.as_deref()),
            check: PssmFileCheck::new(out),
        }
    }

    /// `SavePssmToFile` of `pssm`. An ASN.1 file that cannot be written is LOSAT's refusal there,
    /// where NCBI stops with its serial stream's write error (DECISIONS 10 of line PS: after the
    /// report of the round before; no partial file is left); an ASCII file that cannot be opened
    /// gets nothing, as NCBI's unchecked `ofstream`.
    pub(crate) fn write(&mut self, pssm: &QueryPssm) -> Result<()> {
        if let (Some(stream), Some(bytes)) = (self.asn.as_mut(), pssm.asn.as_deref()) {
            let file = stream.next_file();
            self.check.check("-out_pssm", &file.path)?;
            let written = file.open().and_then(|mut out| {
                out.write_all(bytes)?;
                out.flush()
            });
            if let Err(error) = written {
                if is_regular_file(&file.path) {
                    let _ = std::fs::remove_file(&file.path);
                }
                anyhow::bail!(
                    "a -out_pssm file that cannot be written ({}: {error}), where NCBI BLAST+'s psiblast stops after the round's report with a serial-stream write error, is not supported by LOSAT's PSIBLAST",
                    file.path.display()
                );
            }
        }
        if let (Some(stream), Some(bytes)) = (self.ascii.as_mut(), pssm.ascii.as_deref()) {
            let file = stream.next_file();
            self.check.check("-out_ascii_pssm", &file.path)?;
            if let Ok(mut out) = file.open() {
                let _ = out.write_all(bytes);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set(keys: &[u32]) -> InclusionSet {
        keys.iter().map(|&key| SetKey::Index(key)).collect()
    }

    // psiblast_iteration.cpp:68-96 (quoted at `IterationState`): no convergence before two
    // rounds, an empty set after a non-empty one, and every current id among the previous ones
    // (two empty sets too).
    #[test]
    fn has_converged_follows_ncbis_three_rules() {
        let mut state = IterationState::new(0);
        assert!(!state.has_converged());
        state.advance(set(&[]));
        assert!(!state.has_converged(), "one round: never converged");
        state.advance(set(&[]));
        assert!(state.has_converged(), "two empty sets converge");
        let mut state = IterationState::new(0);
        state.advance(set(&[1, 2]));
        assert!(!state.has_converged());
        state.advance(set(&[]));
        assert!(state.has_converged(), "an empty set after a non-empty one");
        let mut state = IterationState::new(0);
        state.advance(set(&[1, 2]));
        state.advance(set(&[2]));
        assert!(state.has_converged(), "a subset of the previous set");
        let mut state = IterationState::new(0);
        state.advance(set(&[1, 2]));
        state.advance(set(&[2, 3]));
        assert!(!state.has_converged(), "a new id");
        assert!(state.go(), "no limit");
    }

    // psiblast_iteration.cpp:58-65, 113-129: the number of rounds, the round's number and the
    // previous set that the report reads.
    #[test]
    fn the_rounds_run_out_and_number_from_one() {
        let mut state = IterationState::new(2);
        assert!(state.go());
        assert_eq!(state.number(), 1);
        assert!(state.previously_found().is_empty());
        state.advance(set(&[7]));
        assert_eq!(state.number(), 2);
        assert_eq!(state.previously_found(), &set(&[7]));
        assert!(state.go());
        state.advance(set(&[7, 8]));
        assert!(!state.has_more_iterations());
        assert!(!state.go());
        assert!(!state.has_converged());
        let mut state = IterationState::new(2);
        state.advance(set(&[7]));
        state.advance(set(&[7]));
        assert!(state.has_converged(), "converged on the last allowed round");
    }

    #[test]
    fn keys_of_records_and_ids_differ() {
        let ids: InclusionSet = [SetKey::Index(3), SetKey::Id(b"gb|X1.1|".to_vec())]
            .into_iter()
            .collect();
        assert!(ids.contains(&SetKey::Index(3)));
        assert!(!ids.contains(&SetKey::Id(b"3".to_vec())));
    }
}
