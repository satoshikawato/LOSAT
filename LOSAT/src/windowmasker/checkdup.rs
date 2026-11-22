//! `-checkdup`: windowmasker's check of the input for duplicated sequences (`CheckDuplicates`),
//! which `-mk_counts` runs before it counts. Samples of 100 letters every 10,000 letters of each
//! sequence of 50,000 letters or more go into a table; every 100 letters of each later sequence
//! are looked up in it, and runs of at least 4 samples found at the same distances in an earlier
//! sequence are reported as possible duplications (a `LOG_POST` of a `Warning`, without a
//! prefix). The counts do not change.
//!
//! NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_dup_table.cpp:55-586

use rustc_hash::FxHashMap;

use crate::ncbi_app::files::Ostream;
use crate::ncbi_app::Diag;
use crate::seqid::label::both_version;
use crate::seqid::rank::{find_best, text_score};
use crate::seqid::ParsedSeqId;

use super::counts::{post_after_output, GenomeInput};
use super::work::Failure;

/// The length of a sample.
const SAMPLE_LENGTH: usize = 100;
/// The distance between the samples of a sequence.
const SAMPLE_SKIP: u32 = 10_000;
/// The shortest sequence that is checked.
const MIN_SEQ_LENGTH: u32 = 50_000;
/// How far from 10,000 letters the next sample of a run may be found.
const MAX_OFFSET_ERROR: usize = 5;
/// The fewest samples in a run that is reported.
const MIN_MATCH_COUNT: u32 = 4;

/// A place of a sample: the number of its sequence among those checked (from 0, over all
/// files) and its offset (`sample_loc`, ordered by both).
type SampleLoc = (usize, u32);

/// `dup_lookup_table`: the ids of the sequences checked so far and the places of their samples.
/// `next_offset` is the function-static offset of the next sample, which carries from one
/// sequence into the next.
#[derive(Default)]
struct DupTable {
    ids: Vec<Vec<u8>>,
    samples: FxHashMap<Vec<u8>, Vec<SampleLoc>>,
    next_offset: u32,
}

impl DupTable {
    /// `add_seq_info`: the samples of a checked sequence.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_dup_table.cpp:247-266
    fn add_seq_info(&mut self, id: Vec<u8>, data: &[u8]) {
        self.ids.push(id);
        let seqnum = self.ids.len() - 1;
        let data_len = data.len() as u32;
        while self.next_offset < data_len - SAMPLE_LENGTH as u32 {
            let at = self.next_offset as usize;
            self.samples
                .entry(data[at..at + SAMPLE_LENGTH].to_vec())
                .or_default()
                .push((seqnum, self.next_offset));
            self.next_offset += SAMPLE_SKIP;
        }
        self.next_offset = self.next_offset.saturating_sub(data_len);
    }
}

/// A run being followed (`tracker::result`): its number of samples, the place in the earlier
/// sequence where its next sample is expected, and the offset in the current sequence of its
/// last sample.
#[derive(Clone, Copy)]
struct Run {
    count: u32,
    loc: SampleLoc,
    s_offset: usize,
}

/// `tracker`: the runs of the current sequence (`subject_id`) against the earlier ones.
struct Tracker<'a> {
    table: &'a DupTable,
    subject_id: &'a [u8],
    main_list: Vec<Run>,
    aux_list: Vec<Run>,
}

impl Tracker<'_> {
    /// `report_match`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_dup_table.cpp:374-387
    fn report_match(
        &self,
        run: &Run,
        out: &mut Ostream<'_>,
        diag: &mut Diag,
    ) -> Result<(), Failure> {
        let s_off = run.s_offset as u64 + u64::from(SAMPLE_SKIP);
        let q_off = u64::from(run.loc.1);
        let span = u64::from(run.count.wrapping_mul(SAMPLE_SKIP));
        let skip = u64::from(SAMPLE_SKIP);
        let mut text = b"Possible duplication of sequences:\nsubject: ".to_vec();
        text.extend_from_slice(self.subject_id);
        text.extend_from_slice(b" and query: ");
        text.extend_from_slice(&self.table.ids[run.loc.0]);
        text.extend_from_slice(
            format!(
                "\nat intervals\nsubject: {} --- {}\nquery  : {} --- {}\n",
                s_off.wrapping_sub(span),
                s_off.wrapping_sub(skip),
                q_off.wrapping_sub(span),
                q_off.wrapping_sub(skip),
            )
            .as_bytes(),
        );
        post_after_output(out, diag, &text)
    }

    /// The end of a run: reported when it has enough samples.
    fn end(&self, run: &Run, out: &mut Ostream<'_>, diag: &mut Diag) -> Result<(), Failure> {
        if run.count >= MIN_MATCH_COUNT {
            self.report_match(run, out, diag)?;
        }
        Ok(())
    }

    /// A run whose next sample was not at `subject_offset`: ended when the offset has passed
    /// it, else kept.
    fn pass(
        &mut self,
        run: Run,
        subject_offset: usize,
        out: &mut Ostream<'_>,
        diag: &mut Diag,
    ) -> Result<(), Failure> {
        if run.s_offset + SAMPLE_SKIP as usize + MAX_OFFSET_ERROR < subject_offset {
            self.end(&run, out, diag)?;
        } else {
            self.aux_list.push(run);
        }
        Ok(())
    }

    /// `tracker::operator()`: the sample at `subject_offset` was found at `locs` (in order):
    /// the runs and the places are merged in their order; a run whose expected place is found
    /// within 5 letters of 10,000 after its last sample grows, the others end or wait, and each
    /// place starts or continues a run.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_dup_table.cpp:408-500
    fn track(
        &mut self,
        subject_offset: usize,
        locs: &[SampleLoc],
        out: &mut Ostream<'_>,
        diag: &mut Diag,
    ) -> Result<(), Failure> {
        let mut main = std::mem::take(&mut self.main_list);
        let (mut r, mut it) = (0, 0);
        let start = |loc: SampleLoc, count| Run {
            count,
            loc: (loc.0, loc.1.wrapping_add(SAMPLE_SKIP)),
            s_offset: subject_offset,
        };
        loop {
            if r == main.len() {
                if it == locs.len() {
                    break;
                }
                self.aux_list.push(start(locs[it], 1));
                it += 1;
            } else if it == locs.len() || locs[it] > main[r].loc {
                self.pass(main[r], subject_offset, out, diag)?;
                r += 1;
            } else if locs[it] < main[r].loc {
                self.aux_list.push(start(locs[it], 1));
                it += 1;
            } else {
                let mut count = 1;
                while r < main.len() && main[r].loc == locs[it] {
                    let run = main[r];
                    let expected = run.s_offset + SAMPLE_SKIP as usize;
                    if subject_offset < expected - MAX_OFFSET_ERROR {
                        self.aux_list.push(run);
                    } else if subject_offset > expected + MAX_OFFSET_ERROR {
                        self.end(&run, out, diag)?;
                    } else {
                        count = run.count + 1;
                    }
                    r += 1;
                }
                self.aux_list.push(start(locs[it], count));
                it += 1;
            }
        }
        // The old list, emptied, takes the next runs: both lists keep their room, as NCBI's
        // `clear` and `swap` keep it.
        main.clear();
        self.main_list = std::mem::replace(&mut self.aux_list, main);
        Ok(())
    }

    /// `~tracker`: the runs left, reported when they have enough samples.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_dup_table.cpp:390-405
    fn finish(&self, out: &mut Ostream<'_>, diag: &mut Diag) -> Result<(), Failure> {
        for run in &self.main_list {
            self.end(run, out, diag)?;
        }
        Ok(())
    }
}

/// The label of a sequence's best id (`sequence::GetId` with `eGetId_Best`, the first id of
/// the lowest `CSeq_id::Score`, then `GetLabel` in its default form: `lcl|1`,
/// `gnl|BL_ORD_ID:0`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:416-422,498-501
fn best_id_label(ids: &[ParsedSeqId]) -> Vec<u8> {
    let mut label = Vec::new();
    if let Some(id) = find_best(ids, text_score) {
        both_version(id, &mut label);
    }
    label
}

/// `CheckDuplicates` over the input of `-mk_counts`: one read of every file, the reader's
/// posts and the reports in the order NCBI makes them, each after the output written so far.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_dup_table.cpp:540-583
pub(super) fn check_duplicates(input: &mut GenomeInput<'_, '_, '_, '_>) -> Result<(), Failure> {
    let mut table = DupTable::default();
    input.each_record(&mut |record, out, diag| {
        let data = &record.iupac;
        let length = data.len() as u32;
        if length < MIN_SEQ_LENGTH {
            return Ok(());
        }
        let id = best_id_label(&record.ids);
        let data_len = (length - SAMPLE_SKIP) as usize;
        let mut tracker = Tracker {
            table: &table,
            subject_id: &id,
            main_list: Vec::new(),
            aux_list: Vec::new(),
        };
        if !table.samples.is_empty() {
            for i in 0..data_len {
                if let Some(locs) = table.samples.get(&data[i..i + SAMPLE_LENGTH]) {
                    tracker.track(i, locs, out, diag)?;
                }
            }
        }
        let main_list = std::mem::take(&mut tracker.main_list);
        drop(tracker);
        table.add_seq_info(id.clone(), data);
        let tracker = Tracker {
            table: &table,
            subject_id: &id,
            main_list,
            aux_list: Vec::new(),
        };
        tracker.finish(out, diag)
    })
}
