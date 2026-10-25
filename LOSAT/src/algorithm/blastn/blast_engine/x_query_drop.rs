//! EXPERIMENT (LOSAT_X_BNQDROP; in the strict set since 2026-10-25: on unless it is `0`): the
//! letters of a query batch's records dropped once the batch has encoded them, where nothing reads
//! them after the search, instead of holding them until the end of the search.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
//! 	    BLAST_PROF_START( APP.LOOP.PRE );
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//! ```
//! NCBI holds a batch's Bioseqs in its scope until the batch is formatted, beside the search's
//! own encoded copy (`SetupQueries_OMF`). No NCBI counterpart for the drop: it is memory only.
//!
//! What it changes. Where the query records are read before the search (`search_cli` without
//! LOSAT_X_QSTREAM acting), `search_cli` holds every record, with its IUPAC letters, until the
//! search ends, and `search_query_batch` encodes each batch's letters into its BLASTNA buffer (both
//! strands and sentinels), which the scan and the traceback read. After the encoding the letters
//! are read only by the batch's query chunks (`search_query_chunks`: their cut, DUST and encoding)
//! and by the reports that show a query's letters. With the switch, where the rule below holds,
//! `search_cli` moves the records into a `Store` and the batch loop takes each batch's records out
//! of it (`x_query_stream::Queries::Taken`), so the loop holds them alone: `search_query_batch`
//! drops their letters after its last reader (the encoding and the query chunks), before it builds
//! the lookup tables (`record.sequence = Vec::new()`, as LOSAT_X_BNSTORE does with the subjects'
//! letters). The loop keeps the records without letters for a writer after the batch, as the
//! streamed records are kept. A batch split into query chunks cuts them from its records' letters,
//! so it drops the letters after its chunks, before its traceback; the chunks' own copies
//! (`FastaRecord::cut`) keep theirs until the chunk ends, as before (dropping a copy of about 1 MiB
//! before its chunk's scan saved nothing measurable: blastn of the 5.5 Mb query below, whose
//! records LOSAT_X_QSTREAM reads batch by batch, read the same max RSS).
//!
//! The rule (`Rule::acts`; the reach of the switch is recorded where `search_cli` decides): every
//! output format reads no query letter after the search (`BlastnOutputFormat::reads_query_letters`:
//! the tabular formats 6 and 7, the Seq-aligns of 8, 9 and 12, and SAM without both `SQ` and `SR`;
//! not outfmt 0, 5, 11, 13 to 16, nor SAM's `SQ SR`) and no caller takes the hits' records (their
//! query rows), without `-query_loc` (its input records keep their letters), where no
//! LOSAT_X_QBATCHPAR slot may run (a slot searches the records of a batch ahead of the loop), and
//! for the records read before the search only (not a streamed search, and not `run_local`, whose
//! caller holds its records).
//!
//! Why the reports are the same. The records keep their local IDs, titles and messages; every
//! length the search and the reports read after the encoding comes from the lengths taken before
//! it (`query_lengths` in the traceback, the loop's searched lengths for an invalid query's
//! message, `QueryNames::searched_lengths` in the Seq-aligns), and no format of the rule reads a
//! query letter (the HSP rows of outfmt 0 and the XML and JSON formats, the archive's Bioseqs and
//! its DUST, and SAM's SEQ of the query with `SQ SR` do). Dropping memory changes no value that the
//! search computes. No floating-point operation is involved.
//!
//! Memory: the peak of a megablast batch is the scan (with LOSAT_X_BNLUTDROP), where the letters
//! of a 5.5 Mb query (E. coli Sakai x MG1655) were 5.2 MiB beside the 10.5 MiB of the BLASTNA
//! buffer and the 106 MiB of the lookup table (plan v4 line MF, `REPORT.md` §2.3): its peak
//! `RssAnon` falls by 4.4-6.4 MiB at one and 16 threads, `-subject` and `-db`. Megablast of 31,168
//! cDNAs at one thread (37 batches, 89 MB of letters) drops each batch's letters after its
//! encoding; its peak, in the first batches, falls by 4.4 MiB.
//!
//! Timing: with LOSAT_TIMING, `[TIMING] x_bn_q_drop: batches= records= letters=` gives the number
//! of batches and records whose letters were dropped, and their letters.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock, PoisonError};

use crate::algorithm::blastn::hsp::BlastnOutputFormat;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::QueryInput;

/// No NCBI counterpart: LOSAT_X_BNQDROP (read once).
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNQDROP"))
}

/// No NCBI counterpart: what `search_cli` knows of a search when it decides whether the switch
/// acts.
pub(super) struct Rule<'a> {
    pub formats: &'a [BlastnOutputFormat],
    /// `-query_loc` is given.
    pub query_range: bool,
    /// A LOSAT_X_QBATCHPAR slot may search a batch ahead of the loop.
    pub slots_may_run: bool,
    /// The caller takes the hits' records (`ReportOutputs::hits`, which hold the query rows).
    pub hit_records: bool,
}

impl Rule<'_> {
    /// Whether the records read before the search are moved into a `Store` (the switch on).
    pub(super) fn acts(&self) -> bool {
        !self.query_range
            && !self.slots_may_run
            && !self.hit_records
            && !self
                .formats
                .iter()
                .any(|format| format.reads_query_letters())
    }
}

/// The query records read before a search, which the batch loop takes out batch by batch, in
/// order (`x_query_stream::Queries::Taken`).
pub(super) struct Store {
    /// The records the loop has not taken, and how many it has taken.
    left: Mutex<(std::vec::IntoIter<FastaRecord>, usize)>,
    /// Every record searched whole (no `-query_loc`).
    query_input: QueryInput,
}

impl Store {
    pub(super) fn new(records: Vec<FastaRecord>) -> Self {
        Self {
            query_input: QueryInput::whole(&records),
            left: Mutex::new((records.into_iter(), 0)),
        }
    }

    /// The query input (every record searched whole).
    pub(super) fn query_input(&self) -> &QueryInput {
        &self.query_input
    }

    /// The length of every record.
    pub(super) fn lengths(&self) -> &[usize] {
        &self.query_input.input_lengths
    }

    /// The records `range`, taken out: the loop takes its batches in order, each starting where
    /// the one before ended.
    pub(super) fn take(&self, range: std::ops::Range<usize>) -> anyhow::Result<Vec<FastaRecord>> {
        let mut left = self.left.lock().unwrap_or_else(PoisonError::into_inner);
        let (records, taken) = &mut *left;
        anyhow::ensure!(
            range.start == *taken && range.len() <= records.len(),
            "LOSAT_X_BNQDROP: the batch of the records {range:?} after those up to {taken}"
        );
        *taken = range.end;
        Ok(records.by_ref().take(range.len()).collect())
    }

    /// A copy of the records `range` for a LOSAT_X_QBATCHPAR slot; `None` once the loop has taken
    /// the first of them. No slot runs where the loop takes the records (`x_qbatch_par::Off::Taken`;
    /// the rule acts where none may run), so a slot never waits on the loop's current batch.
    pub(super) fn copy(&self, range: std::ops::Range<usize>) -> Option<Vec<FastaRecord>> {
        let left = self.left.lock().unwrap_or_else(PoisonError::into_inner);
        let (records, taken) = &*left;
        let from = range.start.checked_sub(*taken)?;
        Some(records.as_slice()[from..from + range.len()].to_vec())
    }

    /// Writes the reader's messages of every record (none taken yet).
    pub(super) fn write_warnings(
        &self,
        out: &mut (dyn std::io::Write + Send),
    ) -> anyhow::Result<()> {
        let left = self.left.lock().unwrap_or_else(PoisonError::into_inner);
        for record in left.0.as_slice() {
            out.write_all(&record.warnings)?;
        }
        Ok(())
    }
}

/// The query records of a batch as `search_query_batch` reads them: shared with other readers, or
/// held by the caller alone (the loop's batch taken from a `Store`), whose letters the batch drops
/// after its last reader of them.
pub(super) enum BatchQueries<'r> {
    Shared(&'r [FastaRecord]),
    Held(&'r mut [FastaRecord]),
}

impl BatchQueries<'_> {
    pub(super) fn records(&self) -> &[FastaRecord] {
        match self {
            BatchQueries::Shared(records) => records,
            BatchQueries::Held(records) => records,
        }
    }

    /// The batch has encoded its letters and searched its query chunks: the letters of records the
    /// caller holds alone are dropped (with the switch).
    pub(super) fn after_encoding(&mut self) {
        let BatchQueries::Held(records) = self else {
            return;
        };
        if !on() {
            return;
        }
        let mut letters = 0;
        for record in records.iter_mut() {
            letters += record.sequence.len() as u64;
            record.sequence = Vec::new();
        }
        BATCHES.fetch_add(1, Ordering::Relaxed);
        RECORDS.fetch_add(records.len() as u64, Ordering::Relaxed);
        LETTERS.fetch_add(letters, Ordering::Relaxed);
    }
}

// No NCBI counterpart: the batches, records and letters dropped after their encoding.
static BATCHES: AtomicU64 = AtomicU64::new(0);
static RECORDS: AtomicU64 = AtomicU64::new(0);
static LETTERS: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: with LOSAT_TIMING, the `[TIMING] x_bn_q_drop` line `main` prints at exit
/// (when a batch dropped its letters).
pub(super) fn print_summary() {
    let batches = BATCHES.load(Ordering::Relaxed);
    if crate::utils::xtiming::on() && batches > 0 {
        eprintln!(
            "[TIMING] x_bn_q_drop: batches={batches} records={} letters={}",
            RECORDS.load(Ordering::Relaxed),
            LETTERS.load(Ordering::Relaxed)
        );
    }
}
