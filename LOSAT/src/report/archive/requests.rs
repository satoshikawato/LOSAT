//! The request of one query batch (`Blast4-queue-search-request`) as the archive (outfmt 11) and
//! the search strategy file (`-export_search_strategy`) take it, made from the batch's records: one
//! builder for the five programs, with what differs between them (the molecules and the frames of
//! the `LCaseMask` parameters) as the run's `RequestPolicy`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:571-573,618-622
//! `CExportStrategy` makes the request of a batch from its query factory: the batch's Bioseqs and
//! Seq-locs, the `-query_loc` interval of its first query, and the user's masks (the lower case of
//! `-lcase_masking`), with the `-subject` sequences (`x_Process_SearchDb`).

use std::sync::Arc;

use anyhow::Result;

use super::bioseq::BioseqInput;
use super::masks::{Frame, MaskEntry};
use super::strategy::StrategyExport;
use super::writer::{RequestInputs, RequestQueries};
use crate::blastinput::seq_range::RecordPlacement;
use crate::report::serial::pssm::PssmWithParameters;

/// What a request reads of a query or subject record. No NCBI counterpart: LOSAT has two record
/// types (`blastinput::fasta_reader::FastaRecord`, BLASTX's `blastx::input::FastaRecord`).
pub(crate) trait RequestRecord {
    /// The record as the Bioseq of a request, of molecule `protein`, with its defline line (the
    /// line as read, else `>` and the title); a query of `-parse_deflines` with an id that is not
    /// local also as the Seq-loc of its searched interval (`seq_loc`).
    fn bioseq(&self, protein: bool) -> BioseqInput<'_>;
    /// The number of letters searched (those of a `-query_loc` interval).
    fn searched_len(&self) -> usize;
    /// The lower-case runs `[from, to]` of the whole record (those of `-lcase_masking`), in the
    /// record's coordinates.
    fn lowercase_runs(&self) -> Vec<(i32, i32)>;
    /// The mask of the record in frame `frame` with `intervals`.
    fn mask(&self, frame: Frame, intervals: Vec<(i32, i32)>) -> MaskEntry;
}

impl RequestRecord for crate::blastinput::fasta_reader::FastaRecord {
    fn bioseq(&self, protein: bool) -> BioseqInput<'_> {
        BioseqInput {
            id: &self.local_id,
            ids: self.believed().map(|ids| ids.ids.as_slice()),
            seq_loc: self.archive_seq_loc(),
            defline: super::defline_of(self),
            title: &self.title,
            letters: self.archive_letters(),
            protein,
        }
    }

    fn searched_len(&self) -> usize {
        self.sequence.len()
    }

    fn lowercase_runs(&self) -> Vec<(i32, i32)> {
        crate::algorithm::common::masks::lowercase_runs(self.archive_letters())
    }

    fn mask(&self, frame: Frame, intervals: Vec<(i32, i32)>) -> MaskEntry {
        MaskEntry {
            id: self.local_id.clone(),
            best: self.believed().map(|ids| ids.best.clone()),
            frame,
            intervals,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1102-1134
// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1084-1085
/// How the lower-case runs of a query (the user's mask of `-lcase_masking`) become the request's
/// `LCaseMask` parameters: one `Blast4-mask` per query, or per frame change of a translated
/// query, from the user's mask of the query, which holds each run once for a protein query and
/// once for each strand (`plus1`, `minus1`) for a nucleotide query.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum LcaseFrames {
    /// A protein query (BLASTP, TBLASTN): one `notset` mask with every run.
    NotSet,
    /// A nucleotide query that is not translated (BLASTN): one `plus1` mask with every run twice,
    /// for the plus and the minus strand.
    PlusTwice,
    /// A translated query (BLASTX, TBLASTX): the frame changes at every interval of the user's
    /// mask, so each run is a mask in `plus1` and one in `minus1`.
    PlusAndMinus,
}

impl LcaseFrames {
    /// The `LCaseMask` parameters of `records`, in order (a record without lower case has none).
    fn masks<R: RequestRecord>(self, records: &[R]) -> Vec<MaskEntry> {
        let mut masks = Vec::new();
        for record in records {
            let runs = record.lowercase_runs();
            match self {
                LcaseFrames::NotSet if !runs.is_empty() => {
                    masks.push(record.mask(Frame::NotSet, runs));
                }
                LcaseFrames::PlusTwice if !runs.is_empty() => {
                    let twice = runs.iter().flat_map(|&run| [run, run]).collect();
                    masks.push(record.mask(Frame::Plus1, twice));
                }
                LcaseFrames::PlusAndMinus => {
                    for run in runs {
                        for frame in [Frame::Plus1, Frame::Minus1] {
                            masks.push(record.mask(frame, vec![run]));
                        }
                    }
                }
                _ => {}
            }
        }
        masks
    }
}

/// What the requests of a program's run take from the program: the molecules of the queries and
/// of the subjects, and how the lower-case runs become the `LCaseMask` parameters (`None`: the
/// request has none, as without `-lcase_masking`, or the caller writes its own, BLASTN's archive;
/// BLASTX's reader keeps runs only with `-lcase_masking`, so its policy is always `PlusAndMinus`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct RequestPolicy {
    pub query_protein: bool,
    pub subject_protein: bool,
    pub lcase: Option<LcaseFrames>,
}

/// The requests of a run's query batches (or of the batches of one call that writes their
/// archives): the policy, the `-subject` records as the Bioseqs that every request holds (made
/// once; `None` for `-db`) and the PSSM of a PSSM query (`with_pssm_query`).
pub(crate) struct RunRequests<'a> {
    policy: RequestPolicy,
    subjects: Vec<BioseqInput<'a>>,
    /// The PSSM that every request holds in place of the batch's queries (`with_pssm_query`).
    pssm: Option<Arc<PssmWithParameters>>,
}

impl<'a> RunRequests<'a> {
    /// The requests of a run of `policy`, with the `-subject` records `subjects` (`None` for `-db`),
    /// whose Bioseqs are made here once for every request.
    pub(crate) fn new<R: RequestRecord>(policy: RequestPolicy, subjects: Option<&'a [R]>) -> Self {
        RunRequests {
            policy,
            subjects: subjects
                .unwrap_or_default()
                .iter()
                .map(|record| BioseqInput {
                    seq_loc: None,
                    ..record.bioseq(policy.subject_protein)
                })
                .collect(),
            pssm: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:338-340
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:307
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:281-288
    /// The requests of a run whose query is a PSSM (psiblast's and tblastn's `-in_pssm`, `pssm`:
    /// the object as the application read it, `PsiSearch::request`): each request's queries are
    /// the PSSM (`queries pssm`), with no `-query_loc` range (the option excludes `-in_pssm`) and no
    /// `LCaseMask` (`x_Process_Pssm` takes no query masks). `None` keeps the batches' queries.
    pub(crate) fn with_pssm_query(mut self, pssm: Option<Arc<PssmWithParameters>>) -> Self {
        self.pssm = pssm;
        self
    }

    /// Calls `each` with the request of the batch of `queries` (each cut to its `-query_loc`;
    /// `first_placement` is where the first one's letters lie in its record, `None` without a
    /// range): the archive adds its results to it (`BatchArchive`), the search strategy writes it.
    pub(crate) fn with_batch<R: RequestRecord, T>(
        &self,
        queries: &[R],
        first_placement: Option<RecordPlacement>,
        each: impl FnOnce(&RequestInputs<'_>) -> T,
    ) -> T {
        if let Some(pssm) = self.pssm.as_deref() {
            return each(&RequestInputs {
                queries: RequestQueries::Pssm(pssm),
                subjects: &self.subjects,
                required: None,
                lcase_masks: &[],
            });
        }
        let query_inputs: Vec<BioseqInput<'_>> = queries
            .iter()
            .map(|record| record.bioseq(self.policy.query_protein))
            .collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:587-605
        // `RequiredStart` / `RequiredEnd`: the interval of the batch's first query when a range
        // cut it to less than its record.
        let required = queries
            .first()
            .zip(first_placement)
            .and_then(|(query, placement)| {
                let cut = query.searched_len();
                (cut != placement.length)
                    .then_some((placement.offset as i32, (placement.offset + cut) as i32 - 1))
            });
        let lcase_masks = self
            .policy
            .lcase
            .map_or_else(Vec::new, |frames| frames.masks(queries));
        each(&RequestInputs {
            queries: RequestQueries::Bioseqs(&query_inputs),
            subjects: &self.subjects,
            required,
            lcase_masks: &lcase_masks,
        })
    }

    /// Writes the search strategy's request of the batch of `queries` (`with_batch`), whose first
    /// query is the run's query `first_query`, where NCBI writes it (`StrategyExport::write_with`):
    /// a request already written (a search of the batches again) or after the run's stop is not
    /// made. NCBI's posts of the bytes that ASN.1 text replaced in it wait in `strategy` for the
    /// run to write them where the batch starts (`StrategyExport::take_posts`).
    pub(crate) fn write<R: RequestRecord>(
        &self,
        strategy: &StrategyExport,
        first_query: usize,
        queries: &[R],
        first_placement: Option<RecordPlacement>,
    ) -> Result<()> {
        strategy
            .write_with(first_query..first_query + queries.len(), |write| {
                self.with_batch(queries, first_placement, |request| write(request))
            })
            .map(drop)
    }
}

#[cfg(test)]
#[path = "requests_tests.rs"]
mod tests;
