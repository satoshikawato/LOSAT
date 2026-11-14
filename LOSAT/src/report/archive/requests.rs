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

use anyhow::Result;

use super::bioseq::BioseqInput;
use super::masks::{Frame, MaskEntry};
use super::strategy::StrategyExport;
use super::writer::RequestInputs;
use crate::blastinput::seq_range::RecordPlacement;

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
/// archives): the policy and the `-subject` records as the Bioseqs that every request holds
/// (made once; `None` for `-db`).
pub(crate) struct RunRequests<'a> {
    policy: RequestPolicy,
    subjects: Vec<BioseqInput<'a>>,
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
        }
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
        let query_inputs: Vec<BioseqInput<'_>> = queries
            .iter()
            .map(|record| record.bioseq(self.policy.query_protein))
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
mod tests {
    use super::*;
    use crate::blastinput::app::Program;
    use crate::blastinput::fasta_reader::FastaRecord;
    use crate::report::archive::options::OptionInputs;
    use crate::report::archive::OptionsHandle;
    use std::sync::{Arc, Mutex};

    const PROTEIN: RequestPolicy = RequestPolicy {
        query_protein: true,
        subject_protein: true,
        lcase: None,
    };

    #[derive(Clone, Default)]
    struct Sink(Arc<Mutex<Vec<u8>>>);

    impl std::io::Write for Sink {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(bytes);
            Ok(bytes.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    fn export() -> (StrategyExport, Sink) {
        let sink = Sink::default();
        let export = StrategyExport::new(
            Box::new(sink.clone()),
            &OptionInputs::new(Program::Blastp, "blastp"),
            None,
            OptionsHandle::default(),
        );
        (export, sink)
    }

    /// A title that ASN.1 text does not write as it is: the request holds `#` and the byte is
    /// posted (NCBI 2.17.0, probe RG K/k4b/probe case blastp_strat: `Error: (803.7) [blastp]
    /// Blast4-request.body.queue-search.queries.bioseq-set.seq-set.E.seq.descr.E.title`), one
    /// error post that the run's archive saves as a message; a subject's title the same; a
    /// request after the run's stop is not made and posts nothing.
    #[test]
    fn a_title_with_a_replaced_byte_is_written_with_a_post() {
        let good = [FastaRecord::new("Query_1", b"q", b"MKV")];
        let bad = [FastaRecord::new("Query_2", b"b\x7fta", b"MKV")];
        let requests = RunRequests::new::<FastaRecord>(PROTEIN, None);
        let (strategy, sink) = export();
        requests.write(&strategy, 0, &good, None).unwrap();
        strategy.note_searched(0..1, None).unwrap();
        assert!(strategy.take_posts(usize::MAX).is_empty());
        let handle = crate::report::archive::ArchiveHandle::new(
            crate::report::archive::ArchiveEncoding::Text,
            crate::report::archive::ArchiveRequest {
                program: "blastp",
                service: "plain",
                algorithm: Vec::new(),
                database: None,
            },
            OptionsHandle::default(),
        );
        requests.write(&strategy, 1, &bad, None).unwrap();
        strategy.note_searched(1..2, None).unwrap();
        // The posts wait for the batch's start, by its first query.
        assert!(strategy.take_posts(1).is_empty());
        let taken = strategy.take_posts(2);
        assert_eq!(taken.len(), 1);
        assert_eq!(taken[0].0, 1);
        let mut diagnostics = Vec::new();
        taken[0].1.write(&mut diagnostics, Some(&handle)).unwrap();
        let text = String::from_utf8(sink.0.lock().unwrap().clone()).unwrap();
        assert!(text.contains("title \"b#ta\""), "{text}");
        let post = "Error: (803.7) [blastp] Blast4-request.body.queue-search.queries.bioseq-set.seq-set.E.seq.descr.E.title\nBad char [0x7F] in string at line 21\nb\x7fta\n";
        assert_eq!(String::from_utf8(diagnostics.clone()).unwrap(), post);
        assert_eq!(
            handle.take_messages(),
            vec![(
                crate::report::diag_log::ERROR,
                b"Error: (803.7) [blastp] Blast4-request.body.queue-search.queries.bioseq-set.seq-set.E.seq.descr.E.title Bad char [0x7F] in string at line 21 b ta".to_vec()
            )]
        );
        // After a stop no request is made.
        strategy.stop_after(2);
        let written = sink.0.lock().unwrap().len();
        requests.write(&strategy, 2, &bad, None).unwrap();
        assert_eq!(sink.0.lock().unwrap().len(), written);
        assert!(strategy.take_all_posts().is_empty());
        // A subject's title: posted with each request that holds it.
        let subjects = [FastaRecord::new("Subject_1", b"s\x7f", b"MKV")];
        let requests = RunRequests::new(PROTEIN, Some(&subjects[..]));
        let (strategy, sink) = export();
        requests.write(&strategy, 0, &good, None).unwrap();
        let text = String::from_utf8(sink.0.lock().unwrap().clone()).unwrap();
        assert!(text.contains("title \"s#\""), "{text}");
        let mut diagnostics = Vec::new();
        strategy
            .take_all_posts()
            .write(&mut diagnostics, None)
            .unwrap();
        let posted = String::from_utf8(diagnostics).unwrap();
        assert!(
            posted.starts_with(
                "Error: (803.7) [blastp] Blast4-request.body.queue-search.subject.sequences.E.descr.E.title\n"
            ),
            "{posted}"
        );
    }

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
        let requests = RunRequests::new::<FastaRecord>(PROTEIN, None);
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
