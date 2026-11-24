//! psiblast's rounds of a FASTA query (`-num_iterations` other than 1; series PS, item PS5a) and of
//! a PSSM query (`-in_pssm` with `-num_iterations` other than 1 or a PSSM file, item PS5b):
//! the plan of a run (`PsiRoundPlan`), NCBI's iteration state (`IterationState`, a port of
//! `CPsiBlastIterationState`), the inclusion set of a round (`InclusionSet`), what one round's
//! search hands the round driver (`PsiRoundIo`, `RoundOutcome`) and the PSSM files the driver
//! writes between rounds (`RoundFiles`). The driver itself runs blastp's engine once per round per
//! query (`blastp::blast_engine::psi_round`).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-321

use std::collections::BTreeSet;
use std::io::Write as _;

use anyhow::Result;

use super::in_pssm_files::InPssmRounds;
use super::pssm::{
    is_regular_file, OutputFileReset, PsiblastPssmOptions, PssmFileCheck, PssmRequest,
};
use crate::algorithm::blastp::args::BlastpApp;
use crate::algorithm::pssm_search::in_pssm::InPssmFile;

/// The rounds of a psiblast run of FASTA queries or of a PSSM query: `-num_iterations` (0: until
/// convergence), the PSSM engine's options and files (`-pseudocount`, `-out_pssm`,
/// `-out_ascii_pssm`, `-save_each_pssm`), `-save_pssm_after_last_round`, and the input PSSM of
/// `-in_pssm` (`input`).
#[derive(Debug, Clone, PartialEq)]
pub struct PsiRoundPlan {
    pub(crate) num_iterations: u32,
    pub(crate) request: PssmRequest,
    pub(crate) save_last: bool,
    pub(crate) input: Option<InPssmRounds>,
}

impl PsiRoundPlan {
    /// The plan of psiblast's options with rounds (`PsiblastPssmOptions::rounds`; `None`
    /// without), with the PSSM file of `-in_pssm` (`file`) and LOSAT's refusal of a file it cannot
    /// write back as NCBI does, made once the database is opened (`in_pssm_files::rounds_input`).
    pub(crate) fn of(options: &PsiblastPssmOptions, file: Option<InPssmFile>) -> Option<Self> {
        let num_iterations = options.rounds?;
        let request = PssmRequest {
            asn_path: options.out_pssm.clone(),
            ascii_path: options.out_ascii_pssm.clone(),
            save_each: options.save_each_pssm,
            pseudocount: options.pseudocount,
            versions: options.versions.clone(),
        };
        let input = file.map(|file| super::in_pssm_files::rounds_input(file, &request));
        Some(Self {
            num_iterations,
            request,
            save_last: options.save_pssm_after_last_round,
            input,
        })
    }

    /// The PSSM files of the run, for the check before the search (`check_pssm_files`); `None`
    /// without `-out_pssm` and `-out_ascii_pssm`.
    pub(crate) fn file_request(&self) -> Option<&PssmRequest> {
        (self.request.asn_path.is_some() || self.request.ascii_path.is_some())
            .then_some(&self.request)
    }
}

/// The key of a subject in the inclusion sets and the rows of the PSSM's multiple alignment
/// (`CSeq_id_Handle` of its Seq-align id, DESIGN_ps5 §2): the index of a `-subject` record
/// (`lcl|Subject_N`), the OID of a database sequence without parsed Seq-ids
/// (`gnl|BL_ORD_ID|oid`, `seqdb::Oid`; not its rank in a restricted search), or a parsed Seq-align id (copies of one record share it). Two parsed ids
/// are one key where they are equal; others that NCBI's handles can take for one are refused
/// (`check_handles`), so the keys are the handles of the runs LOSAT accepts.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum SetKey {
    Index(u32),
    Id(crate::seqid::model::ParsedSeqId),
}

impl SetKey {
    /// The Seq-align id of an index key (`gnl|BL_ORD_ID|oid`, as NCBI names a database sequence
    /// without parsed Seq-ids; a `-subject` record's key is never compared with a parsed id: the
    /// keys of `-subject` records are their indices, with `-parse_deflines` too, whose reader
    /// refuses the ids that NCBI's handles would take for one; RG DECISIONS 27).
    fn ordinal_id(index: u32) -> crate::seqid::model::ParsedSeqId {
        use crate::seqid::model::{Dbtag, ObjectId, ParsedSeqId};
        ParsedSeqId::General(Dbtag {
            db: b"BL_ORD_ID".to_vec(),
            tag: ObjectId::Id(index as i32),
        })
    }

    /// `CSeq_id::Match` of the two keys' Seq-align ids (`Compare` gives `e_YES`).
    pub(crate) fn matches(&self, other: &SetKey) -> bool {
        use crate::seqid::defline_set::{compare, Sic};
        match (self, other) {
            (SetKey::Index(a), SetKey::Index(b)) => a == b,
            (SetKey::Id(a), SetKey::Id(b)) => compare(a, b) == Sic::Yes,
            (SetKey::Index(index), SetKey::Id(id)) | (SetKey::Id(id), SetKey::Index(index)) => {
                compare(&Self::ordinal_id(*index), id) == Sic::Yes
            }
        }
    }

    /// The key's Seq-align id as FASTA, for LOSAT's messages (an ordinal key named as such).
    fn label(&self) -> String {
        let mut out = Vec::new();
        match self {
            SetKey::Index(index) => {
                crate::seqid::label::fasta(&Self::ordinal_id(*index), &mut out);
                out.extend_from_slice(b" (a sequence without parsed Seq-ids)");
            }
            SetKey::Id(id) => crate::seqid::label::fasta(id, &mut out),
        }
        String::from_utf8_lossy(&out).into_owned()
    }

    /// LOSAT's refusal of the keys `a` and `b`, which differ but which NCBI's handles can take for
    /// one sequence (`HandleIndex`), in a run of `app` (psiblast, or deltablast's rounds, whose
    /// inclusion sets and PSSMs are psiblast's).
    pub(crate) fn refuse_pair(a: &SetKey, b: &SetKey, app: BlastpApp) -> anyhow::Error {
        anyhow::anyhow!(
            "the Seq-align ids {} and {} of two database sequences that {} compares (in its inclusion sets, the rows of its PSSM or its two-part description table), which NCBI BLAST+'s Seq-id handles can take for one sequence (ids that differ only in letter case, one accession under gb, emb and dbj, and other ids that LOSAT does not tell apart as NCBI's handles do), are not supported by LOSAT's {}",
            a.label(),
            b.label(),
            app.ncbi_name(),
            app.name()
        )
    }

    /// `handle_superset_key` of the key's Seq-align id.
    fn superset(&self) -> Vec<u8> {
        let mut out = Vec::new();
        match self {
            SetKey::Index(index) => {
                crate::seqid::defline_set::handle_superset_key(&Self::ordinal_id(*index), &mut out)
            }
            SetKey::Id(id) => crate::seqid::defline_set::handle_superset_key(id, &mut out),
        }
        out
    }

    /// LOSAT's refusal of two keys of `keys` that differ but that NCBI's handles can take for
    /// one sequence, and of a `giim` key (`HandleIndex::new`), in a run of `app`.
    pub(crate) fn check_handles<'a>(
        keys: impl IntoIterator<Item = &'a SetKey>,
        app: BlastpApp,
    ) -> Result<()> {
        HandleIndex::new(keys, app).map(|_| ())
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:2926-2938
// A giim id is never one handle with an equal id (`FindOrCreate` makes a new one each time), so
// NCBI's inclusion sets never find it again, every included HSP of it counts as a row, and its
// scope never finds its Bioseq (review 1 of line DBID pssmrows, H1): LOSAT refuses a giim key
// that psiblast compares.
/// The keys that one round of psiblast compares by NCBI's handles (`CSeq_id_Handle`), by their
/// `seqid::defline_set::handle_superset_key`: building it refuses two keys that differ but share
/// that key (NCBI's handles can be one or two) and a `giim` key (built for the rows of a PSSM, and
/// for a previous set, which was such rows); `unsure` then names a key of the index that another
/// key differs from but can be one handle with. `app`: the run's program (psiblast, or deltablast),
/// which the refusals name.
pub(crate) struct HandleIndex<'a> {
    by_superset: std::collections::BTreeMap<Vec<u8>, &'a SetKey>,
}

impl<'a> HandleIndex<'a> {
    pub(crate) fn new(keys: impl IntoIterator<Item = &'a SetKey>, app: BlastpApp) -> Result<Self> {
        let mut by_superset = std::collections::BTreeMap::new();
        for key in keys {
            if let SetKey::Id(crate::seqid::model::ParsedSeqId::Giim(_)) = key {
                anyhow::bail!(
                    "the Seq-align id {} of a database sequence in {}'s PSSM, a giim id, which NCBI BLAST+'s Seq-id handles never take for one sequence with an equal id (it counts each of its alignments as a row of the PSSM's multiple alignment and cannot read its letters), is not supported by LOSAT's {}",
                    key.label(),
                    app.ncbi_name(),
                    app.name()
                );
            }
            match by_superset.entry(key.superset()) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(key);
                }
                std::collections::btree_map::Entry::Occupied(entry) if *entry.get() != key => {
                    return Err(SetKey::refuse_pair(entry.get(), key, app));
                }
                std::collections::btree_map::Entry::Occupied(_) => {}
            }
        }
        Ok(Self { by_superset })
    }

    /// The key of the index that `key` differs from but that NCBI's handles can take for one
    /// with it; `None` where `key` is in the index or shares its superset key with none.
    pub(crate) fn unsure(&self, key: &SetKey) -> Option<&'a SetKey> {
        self.by_superset
            .get(&key.superset())
            .copied()
            .filter(|indexed| *indexed != key)
    }
}

/// A round's inclusion set (`CPsiBlastIterationState::TSeqIds`).
pub(crate) type InclusionSet = BTreeSet<SetKey>;

/// `CPsiBlastIterationState`.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_iteration.cpp:49-131
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
    /// The Seq-id of the next PSSM's query: the PSSM query's of `-in_pssm`, `None` for a FASTA
    /// query (`Query_<n>`).
    pub(crate) query_id: Option<crate::seqid::ParsedSeqId>,
    /// deltablast `-num_iterations 1`: the round writes the report of a search without rounds
    /// (no "Results from round" line, no `# Iteration:` line), whose PSSM files and next PSSM
    /// the round driver still makes.
    pub(crate) silent: bool,
    /// deltablast's round 1, a result after the batch's first (`CDeltaBlast::Run` makes the
    /// results of the batch's queries one after another, and `CPsiBl2Seq`'s search of `-subject`
    /// one for each subject; the round's inclusion set and next PSSM are made of the first,
    /// `(*results)[0]`): this call reports its query and subject only.
    pub(crate) extra_result: bool,
    /// deltablast: the letters of the next PSSM's query, the query of the domains' PSSM (the
    /// record's with the lower case of `-lcase_masking` as X, `delta_hooks::CddPssm`); `None` for
    /// psiblast (the record's).
    pub(crate) pssm_letters: Option<Vec<u8>>,
    /// deltablast from round 2 on: round 1's ancillary data (of the first query of the batch),
    /// with which NCBI makes the last PSSM of `-save_pssm_after_last_round`
    /// (deltablast_app.cpp:569-577 reads the function's round-1 result, `results_1st_query`).
    pub(crate) first_round_ancillary: Option<super::pssm::PssmAncillary>,
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
/// (`GetSeqIds`), the PSSM made from its Seq-aligns where NCBI makes one (after `Advance`,
/// when `itr` continues or `save_last` is set), with NCBI's error of the PSSM engine, which comes after the round's
/// report, and for a round that searched a PSSM the blocks that NCBI writes back into it
/// (`kbp_psi[0]`, `kbp_gap_psi[0]`; `in_pssm_files::InPssmRounds::written_back`).
pub(crate) struct RoundOutcome {
    pub(crate) has_alignments: bool,
    pub(crate) set: InclusionSet,
    pub(crate) next: Option<Result<super::pssm::QueryPssm>>,
    /// NCBI's messages of the next PSSM's reads (`pssm::PssmReadMessages`: `Failed to retrieve
    /// sequence` for each Seq-align whose read fails, the taxonomy texts of its loads), which
    /// come after the round's report, before the PSSM's error.
    pub(crate) next_warnings: super::pssm::PssmReadMessages,
    pub(crate) searched: Option<(crate::core::psi::KarlinBlk, crate::core::psi::KarlinBlk)>,
    /// deltablast: the round's ancillary data (`PsiRoundIo::first_round_ancillary` of the later
    /// rounds, and the data that the ASCII PSSM after the last PSSM prints).
    pub(crate) ancillary: Option<super::pssm::PssmAncillary>,
    /// deltablast from round 2 on, where another round follows and the last PSSM is saved: the
    /// last PSSM that NCBI saves if the next round has no alignments (this round's Seq-aligns
    /// with round 1's ancillary data), with its messages.
    pub(crate) last: Option<(
        Result<super::pssm::QueryPssm>,
        super::pssm::PssmReadMessages,
    )>,
}

/// The PSSM files of a run (`SavePssmToFile`), written between the rounds: each option's
/// `CAutoOutputFileReset` (one run-level `F.N` counter each), the ASN.1 file then the ASCII one.
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:385-407 (cited at
// `pssm::PssmFiles`)
pub(crate) struct RoundFiles {
    asn: Option<OutputFileReset>,
    ascii: Option<OutputFileReset>,
    check: PssmFileCheck,
}

impl RoundFiles {
    /// `app`: the program's name in the refusals (`BlastpApp::ncbi_name`).
    pub(crate) fn new(
        request: &PssmRequest,
        out: Option<&std::path::Path>,
        app: &'static str,
    ) -> Self {
        let [asn, ascii] = request.streams();
        Self {
            asn,
            ascii,
            check: PssmFileCheck::new(out, app),
        }
    }

    /// `SavePssmToFile` of a PSSM whose ASN.1 text is `asn` and ASCII PSSM `ascii` (each as
    /// requested). An ASN.1 file that cannot be written is LOSAT's refusal there, where NCBI stops
    /// with its serial stream's write error (DECISIONS 10 of line PS: after the report of the round
    /// before, or with `-in_pssm` before round 1, after the prolog; no partial file is left); an ASCII file that cannot be opened gets nothing, as
    /// NCBI's unchecked `ofstream`.
    pub(crate) fn write(&mut self, asn: Option<&[u8]>, ascii: Option<&[u8]>) -> Result<()> {
        if let (Some(stream), Some(bytes)) = (self.asn.as_mut(), asn) {
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
                if self.check.app() != "psiblast" {
                    anyhow::bail!(
                        "a -out_pssm file that cannot be written ({}: {error}) is not supported by LOSAT's {}",
                        file.path.display(),
                        self.check.app().to_ascii_uppercase()
                    );
                }
                anyhow::bail!(
                    "a -out_pssm file that cannot be written ({}: {error}), where NCBI BLAST+'s psiblast stops with a serial-stream write error (after the report of the round before; with -in_pssm or -in_msa before round 1), is not supported by LOSAT's PSIBLAST",
                    file.path.display()
                );
            }
        }
        if let (Some(stream), Some(bytes)) = (self.ascii.as_mut(), ascii) {
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
#[path = "rounds_tests.rs"]
mod tests;
