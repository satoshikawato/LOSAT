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

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use anyhow::Result;

use super::bioseq::BioseqInput;
use super::masks::{Frame, MaskEntry};
use super::strategy::StrategyExport;
use super::writer::{RequestInputs, RequestQueries};
use super::ArchiveEncoding;
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
    /// The title, which the request writes as a VisibleString.
    fn title(&self) -> &[u8];
    /// The record's local id (`Query_N`, `Subject_N`), which a refusal names.
    fn local_id(&self) -> &str;
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
/// once; `None` for `-db`), the PSSM of a PSSM query (`with_pssm_query`) and the program's name in
/// a refusal.
pub(crate) struct RunRequests<'a, R> {
    policy: RequestPolicy,
    label: &'static str,
    subject_records: &'a [R],
    subjects: Vec<BioseqInput<'a>>,
    /// Whether the search strategy's first request checked the subjects' titles (`write`).
    subjects_checked: AtomicBool,
    /// The PSSM that every request holds in place of the batch's queries (`with_pssm_query`).
    pssm: Option<Arc<PssmWithParameters>>,
}

impl<'a, R: RequestRecord> RunRequests<'a, R> {
    /// The requests of a run of `policy`, with the `-subject` records `subjects` (`None` for `-db`),
    /// whose Bioseqs are made here once for every request.
    pub(crate) fn new(
        policy: RequestPolicy,
        subjects: Option<&'a [R]>,
        label: &'static str,
    ) -> Self {
        let subject_records = subjects.unwrap_or_default();
        RunRequests {
            policy,
            label,
            subject_records,
            subjects: subject_records
                .iter()
                .map(|record| BioseqInput {
                    seq_loc: None,
                    ..record.bioseq(policy.subject_protein)
                })
                .collect(),
            subjects_checked: AtomicBool::new(false),
            pssm: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:338-340
    // ```c++
    //     CRef<IQueryFactory> query_factory(s_MakeQueryFactory(pssm, query));
    //
    //     SaveSearchStrategy(args, m_CmdLineArgs, query_factory, opts_hndl, pssm, kNumIterations);
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:307
    // ```c++
    //             SaveSearchStrategy(args, m_CmdLineArgs, query_factory, m_OptsHndl, pssm);
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:281-287
    // ```c++
    //         if (pssm.NotEmpty())
    //         {
    //             _ASSERT(queries.Empty());
    //             if(num_iters != 0)
    //             	retval.Reset(new blast::CExportStrategy(pssm, opts_hndl, search_db, client_id, num_iters));
    //             else
    //             	retval.Reset(new blast::CExportStrategy(pssm, opts_hndl, search_db, client_id));
    // ```
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
    pub(crate) fn with_batch<T>(
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
            queries: RequestQueries::Bioseqs(&query_inputs),
            subjects: &self.subjects,
            required,
            lcase_masks: &lcase_masks,
        })
    }

    /// Refuses the titles of `queries` and, the first time, of the subjects that the search
    /// strategy's requests hold when ASN.1 text cannot write them as NCBI does (`check_titles_in`).
    /// A caller whose batch's request is certain to follow output it is about to write (the
    /// reports of the batches before) calls it before that output.
    pub(crate) fn check_titles(&self, queries: &[R]) -> Result<()> {
        const OPTION: &str = "-export_search_strategy";
        super::check_titles_in(ArchiveEncoding::Text, queries, "query", self.label, OPTION)?;
        if !self.subjects_checked.load(Ordering::Relaxed) {
            super::check_titles_in(
                ArchiveEncoding::Text,
                self.subject_records,
                "subject",
                self.label,
                OPTION,
            )?;
            self.subjects_checked.store(true, Ordering::Relaxed);
        }
        Ok(())
    }

    /// Writes the search strategy's request of the batch of `queries` (`with_batch`), whose first
    /// query is the run's query `first_query`, where NCBI writes it (`StrategyExport::write_with`):
    /// a request already written (a search of the batches again) or after the run's stop is not
    /// made, and its titles are not looked at. The titles the request holds are refused first
    /// (`check_titles`); a refused request writes nothing.
    pub(crate) fn write(
        &self,
        strategy: &StrategyExport,
        first_query: usize,
        queries: &[R],
        first_placement: Option<RecordPlacement>,
    ) -> Result<()> {
        strategy
            .write_with(first_query..first_query + queries.len(), |write| {
                self.check_titles(queries)?;
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
            "BLASTP",
            OptionsHandle::default(),
        )
        .unwrap();
        (export, sink)
    }

    /// A query title that ASN.1 text cannot write refuses its request: nothing is written and the
    /// run's cursor stays where it was; a request after the run's stop is not made, so its titles
    /// are not looked at; a subject's title is refused with the first request.
    #[test]
    fn a_refused_title_writes_nothing_and_a_stopped_request_is_not_checked() {
        let good = [FastaRecord::new("Query_1", b"q", b"MKV")];
        let bad = [FastaRecord::new("Query_2", b"b\x7fta", b"MKV")];
        let requests = RunRequests::<FastaRecord>::new(PROTEIN, None, "BLASTP");
        let (strategy, sink) = export();
        requests.write(&strategy, 0, &good, None).unwrap();
        strategy.note_searched(0..1, None).unwrap();
        let written = sink.0.lock().unwrap().len();
        assert!(written > 0);
        let refusal = requests.write(&strategy, 1, &bad, None).unwrap_err();
        assert!(refusal.is::<super::super::TitleRefusal>(), "{refusal}");
        assert!(refusal.to_string().contains("Query_2"), "{refusal}");
        assert_eq!(sink.0.lock().unwrap().len(), written);
        assert!(strategy.takes(1..2));
        // After a stop the bad title is never written, so it is not refused.
        strategy.stop_after(1);
        requests.write(&strategy, 1, &bad, None).unwrap();
        assert_eq!(sink.0.lock().unwrap().len(), written);
        // A subject title: refused with the first request.
        let subjects = [FastaRecord::new("Subject_1", b"s\x7f", b"MKV")];
        let requests = RunRequests::new(PROTEIN, Some(&subjects[..]), "BLASTP");
        let (strategy, sink) = export();
        let refusal = requests.write(&strategy, 0, &good, None).unwrap_err();
        assert!(
            refusal.to_string().contains("subject Subject_1"),
            "{refusal}"
        );
        assert!(sink.0.lock().unwrap().is_empty());
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

    /// A PSSM query's request holds the PSSM in place of the batch's queries, without a range or
    /// `LCaseMask` parameters (`x_Process_Pssm`), whatever the policy and the records.
    #[test]
    fn a_pssm_query_replaces_the_queries_without_range_or_masks() {
        let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(
            &std::fs::read(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/fixtures/db_search/db.ps.asn.rnd.6.asn"
            ))
            .unwrap(),
        )
        .unwrap();
        let requests = RunRequests::<FastaRecord>::new(
            RequestPolicy {
                lcase: Some(LcaseFrames::NotSet),
                ..PROTEIN
            },
            None,
            "PSIBLAST",
        )
        .with_pssm_query(Some(std::sync::Arc::new(pssm.clone())));
        let queries = [FastaRecord::new("Query_1", b"q", b"mkvMKV")];
        let cut = RecordPlacement {
            offset: 4,
            length: 10,
        };
        requests.with_batch(&queries, Some(cut), |request| {
            assert!(matches!(request.queries, RequestQueries::Pssm(held) if *held == pssm));
            assert_eq!(request.required, None);
            assert!(request.lcase_masks.is_empty());
        });
    }

    #[test]
    fn a_cut_first_query_gives_the_required_interval() {
        let requests = RunRequests::<FastaRecord>::new(PROTEIN, None, "BLASTP");
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
