//! The request of one query batch (`Blast4-queue-search-request`) as the archive (outfmt 11) and
//! the search strategy file (`-export_search_strategy`) take it, made from the batch's records: one
//! builder for the five programs, with what differs between them (the molecules and the frames of
//! the `LCaseMask` parameters) as the run's `RequestPolicy`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:571-573,618-622
//! ```c++
//!     CRef<IRemoteQueryData> remote_query(query->MakeRemoteQueryData());
//!     CRef<CBioseq_set> bioseq_set = remote_query->GetBioseqSet();
//!     IRemoteQueryData::TSeqLocs seqloc_list = remote_query->GetSeqLocs();
//!     ...
//!     CObjMgr_QueryFactory* objmgrqf = dynamic_cast<CObjMgr_QueryFactory*>(&*query);
//!     if ( NULL != objmgrqf )
//!     {
//!     	TSeqLocInfoVector user_specified_masks = objmgrqf->ExtractUserSpecifiedMasks();
//! ```
//! `CExportStrategy` makes the request of a batch from its query factory: the batch's Bioseqs and
//! Seq-locs, the `-query_loc` interval of its first query, and the user's masks (the lower case of
//! `-lcase_masking`), with the `-subject` sequences (`x_Process_SearchDb`).

use std::borrow::Cow;

use anyhow::Result;

use super::bioseq::BioseqInput;
use super::masks::{Frame, MaskEntry};
use super::strategy::StrategyExport;
use super::writer::RequestInputs;
use super::ArchiveEncoding;
use crate::blastinput::seq_range::RecordPlacement;

/// What a request reads of a query or subject record. No NCBI counterpart: LOSAT has two record
/// types (`blastinput::fasta_reader::FastaRecord`, BLASTX's `blastx::input::FastaRecord`).
pub(crate) trait RequestRecord {
    /// The record as the Bioseq of a request, whose defline line is `defline` (`defline`) and
    /// molecule `protein`; a query of `-parse_deflines` with an id that is not local also as the
    /// Seq-loc of its searched interval (`seq_loc`).
    fn bioseq<'a>(&'a self, defline: &'a [u8], protein: bool) -> BioseqInput<'a>;
    /// The defline line of the Bioseq: the line as read, else `>` and the title.
    fn defline(&self) -> Cow<'_, [u8]>;
    /// The number of letters searched (those of a `-query_loc` interval).
    fn searched_len(&self) -> usize;
    /// The lower-case runs `[from, to]` of the whole record (those of `-lcase_masking`), in the
    /// record's coordinates.
    fn lowercase_runs(&self) -> Vec<(i32, i32)>;
    /// The mask of the record in frame `frame` with `intervals`.
    fn mask(&self, frame: Frame, intervals: Vec<(i32, i32)>) -> MaskEntry;
    /// The title, which the request writes as a VisibleString.
    fn title(&self) -> &[u8];
    /// The record's local id (`Query_N`, `Subject_N`), which a refusal names.
    fn local_id(&self) -> &str;
}

impl RequestRecord for crate::blastinput::fasta_reader::FastaRecord {
    fn bioseq<'a>(&'a self, defline: &'a [u8], protein: bool) -> BioseqInput<'a> {
        BioseqInput {
            id: &self.local_id,
            ids: self.believed().map(|ids| ids.ids.as_slice()),
            seq_loc: self.archive_seq_loc(),
            defline,
            title: &self.title,
            letters: self.archive_letters(),
            protein,
        }
    }

    fn defline(&self) -> Cow<'_, [u8]> {
        super::defline_of(self)
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

    fn title(&self) -> &[u8] {
        &self.title
    }

    fn local_id(&self) -> &str {
        &self.local_id
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1102-1134
// ```c++
//     ITERATE(TSeqLocInfoVector, query_masks, masking_locations) {
//         CRef<CPacked_seqint> packed_seqint(new CPacked_seqint);
//         ...
//         int current_frame = query_masks->front()->GetFrame();
//         ITERATE(TMaskedQueryRegions, mask_locs, *query_masks) {
//               if  (Blast_QueryIsTranslated(program) && current_frame != (*mask_locs)->GetFrame())
//               {
//                   if (!packed_seqint.Empty())
//                   {
//                      CRef<CBlast4_mask> network_mask = s_CreateBlastMask(*packed_seqint, program);
//                      network_mask->SetFrame(FrameNumber2NetworkFrame(current_frame, program));
//                      retval.push_back(network_mask);
//                   }
//                   current_frame = (*mask_locs)->GetFrame();
//                   packed_seqint.Reset(new CPacked_seqint);
//               }
//
//               packed_seqint->AddInterval((*mask_locs)->GetSeqId(),
//                              (*mask_locs)->GetInterval().GetFrom(),
//                              (*mask_locs)->GetInterval().GetTo());
//         }
//
//         if (!packed_seqint.Empty())
//         {
//              CRef<CBlast4_mask> network_mask = s_CreateBlastMask(*packed_seqint, program);
//              if (Blast_QueryIsTranslated(program))
//                   network_mask->SetFrame(FrameNumber2NetworkFrame(current_frame, program));
//              retval.push_back(network_mask);
//         }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1084-1085
// ```c++
//     /// The frame can only be notset for protein queries or plus1 for
//     /// nucleotide queries
// ```
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
/// of the subjects, and the `LCaseMask` parameters (`None`: none, without `-lcase_masking`, or
/// written by the caller, BLASTN's archive).
#[derive(Clone, Copy, Debug)]
pub(crate) struct RequestPolicy {
    pub query_protein: bool,
    pub subject_protein: bool,
    pub lcase: Option<LcaseFrames>,
}

/// The requests of a run's query batches: the policy, the `-subject` records (every request holds
/// each of them whole; `None` for `-db`) and the program's name in a refusal.
pub(crate) struct RunRequests<'a, R> {
    pub policy: RequestPolicy,
    pub subjects: Option<&'a [R]>,
    pub label: &'static str,
}

impl<R: RequestRecord> RunRequests<'_, R> {
    /// Calls `each` with the request of the batch of `queries` (each cut to its `-query_loc`;
    /// `first_placement` is where the first one's letters lie in its record, `None` without a
    /// range): the archive adds its results to it (`BatchArchive`), the search strategy writes it.
    pub(crate) fn with_batch<T>(
        &self,
        queries: &[R],
        first_placement: Option<RecordPlacement>,
        each: impl FnOnce(&RequestInputs<'_>) -> T,
    ) -> T {
        // The Bioseq inputs borrow their deflines: keep them alive for the whole call.
        let subject_records = self.subjects.unwrap_or_default();
        let subject_deflines: Vec<Cow<'_, [u8]>> = subject_records.iter().map(R::defline).collect();
        let subjects: Vec<BioseqInput<'_>> = subject_records
            .iter()
            .zip(&subject_deflines)
            .map(|(record, defline)| BioseqInput {
                seq_loc: None,
                ..record.bioseq(defline, self.policy.subject_protein)
            })
            .collect();
        let deflines: Vec<Cow<'_, [u8]>> = queries.iter().map(R::defline).collect();
        let query_inputs: Vec<BioseqInput<'_>> = queries
            .iter()
            .zip(&deflines)
            .map(|(record, defline)| record.bioseq(defline, self.policy.query_protein))
            .collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:587-605
        // ```c++
        //         if (seqloc_list.front()->IsInt())
        //         {
        //             const int kStart((int)seqloc_list.front()->GetStart(eExtreme_Positional));
        //             const int kStop((int)seqloc_list.front()->GetStop(eExtreme_Positional));
        //             const int kRangeLength = kStop - kStart + 1;
        //             ...
        //             const int kFullLength =
        //                 bioseq_set->GetSeq_set().front()->GetSeq().GetInst().GetLength();
        //
        //             if (kFullLength != kRangeLength)
        //             {
        //             	x_AddParameterToProgramOptions(CBlast4Field::Get(eBlastOpt_RequiredStart), kStart);
        //             	x_AddParameterToProgramOptions(CBlast4Field::Get(eBlastOpt_RequiredEnd), kStop);
        // ```
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
            queries: &query_inputs,
            subjects: &subjects,
            required,
            lcase_masks: &lcase_masks,
        })
    }

    /// Writes the search strategy's request of the batch of `queries` (`with_batch`), whose first
    /// query is the run's query `first_query`, where NCBI writes it (`StrategyExport::write_request`):
    /// a request already written (a search of the batches again) or after the run's stop is not
    /// made. The titles the request holds are refused first when ASN.1 text cannot write them as
    /// NCBI does (`check_titles_in`).
    pub(crate) fn write(
        &self,
        strategy: &StrategyExport,
        first_query: usize,
        queries: &[R],
        first_placement: Option<RecordPlacement>,
    ) -> Result<()> {
        let span = first_query..first_query + queries.len();
        if !strategy.takes(span.clone()) {
            return Ok(());
        }
        const OPTION: &str = "-export_search_strategy";
        super::check_titles_in(ArchiveEncoding::Text, queries, "query", self.label, OPTION)?;
        if let Some(subjects) = self.subjects {
            super::check_titles_in(
                ArchiveEncoding::Text,
                subjects,
                "subject",
                self.label,
                OPTION,
            )?;
        }
        self.with_batch(queries, first_placement, |request| {
            strategy.write_request(span, request)
        })
        .map(drop)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::fasta_reader::FastaRecord;

    #[test]
    fn the_lower_case_runs_take_the_frames_of_the_program() {
        let records = [
            FastaRecord::new("Query_1", b"a", b"acGTtt"),
            FastaRecord::new("Query_2", b"b", b"ACGT"),
        ];
        let frames = |lcase: LcaseFrames| -> Vec<(Frame, Vec<(i32, i32)>)> {
            lcase
                .masks(&records)
                .into_iter()
                .map(|mask| {
                    assert_eq!(mask.id, "Query_1");
                    (mask.frame, mask.intervals)
                })
                .collect()
        };
        assert_eq!(
            frames(LcaseFrames::NotSet),
            [(Frame::NotSet, vec![(0, 1), (4, 5)])]
        );
        assert_eq!(
            frames(LcaseFrames::PlusTwice),
            [(Frame::Plus1, vec![(0, 1), (0, 1), (4, 5), (4, 5)])]
        );
        assert_eq!(
            frames(LcaseFrames::PlusAndMinus),
            [
                (Frame::Plus1, vec![(0, 1)]),
                (Frame::Minus1, vec![(0, 1)]),
                (Frame::Plus1, vec![(4, 5)]),
                (Frame::Minus1, vec![(4, 5)]),
            ]
        );
    }

    #[test]
    fn a_cut_first_query_gives_the_required_interval() {
        let requests = RunRequests::<FastaRecord> {
            policy: RequestPolicy {
                query_protein: true,
                subject_protein: true,
                lcase: None,
            },
            subjects: None,
            label: "BLASTP",
        };
        let queries = [FastaRecord::new("Query_1", b"q", b"MKV")];
        let required = |placement| requests.with_batch(&queries, placement, |r| r.required);
        assert_eq!(required(None), None);
        let whole = RecordPlacement {
            offset: 0,
            length: 3,
        };
        assert_eq!(required(Some(whole)), None);
        let cut = RecordPlacement {
            offset: 4,
            length: 10,
        };
        assert_eq!(required(Some(cut)), Some((4, 6)));
    }
}
