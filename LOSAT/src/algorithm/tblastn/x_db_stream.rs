//! EXPERIMENT (LOSAT_X_TBNDBSTREAM / LOSAT_X_TBNDBSTREAMSHADOW): a `-db` search of TBLASTN reads
//! each database sequence where it searches it, instead of a copy of the whole database made
//! before the first query batch.
//!
//! Without the switch, `TblastnDatabase::read_stores` packs every sequence of the database
//! (`common::database::read_nucleotide_stores`) into anonymous memory before the first query
//! batch, and every stage reads the copy until the reports: 23.9 MiB for the 100 Mb of the
//! C. elegans genome, 32.2 MiB for a zebrafish cDNA set, alive at the peak of every thread count
//! (v4-D12 s8 probe, workstation), where NCBI maps the sequence file. The reports then cloned the
//! packed sequences they show (`database::reported_subjects`).
//!
//! With the switch on:
//!
//! - `XDbStream::new` reads the database once in the blocks of `read_nucleotide_stores`, with its
//!   checks in its order (`common::database::x_check_nucleotide_stores`), and keeps the length of
//!   every sequence (four bytes each). Its error is the error of `read_stores`, which fails the
//!   same query batch with the same text.
//! - The preliminary stage of a subject reads its sequence when it makes the subject's plan
//!   (`search_seed::x_scan_plan`, `XPrelimSubject::Stream`) and drops it once the plan is made:
//!   the frames are translated from it and the stage reads nothing else of the sequence. A plan
//!   that the search keeps between query batches (`XSubjectSide`) is not made again, and its
//!   sequence is not read again.
//! - The composition redo and the traceback of `-comp_based_stats 0` read the sequence of each
//!   subject they align (`SubjectLetters::x_nt`), and the reports read the sequences they show
//!   (`database::reported_subjects`), one read each, as NCBI fetches the sequences of its hit
//!   lists.
//!
//! Every read packs the stored bytes as the copy does (`common::database::x_read_nucleotide_store`,
//! the same `PackedNucleotide::from_stored`), so every stage computes the same values on the same
//! subjects in the same order. A read that fails during the search (a database changed or damaged
//! after the check) is the error of `read_nucleotide_stores` for the same sequence.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:375 (`GetSequence` points
//! into the mapped volume) and c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:328-350 (the
//! threads read the database in chunks of consecutive ordinal ids).
//!
//! `LOSAT_X_TBNDBSTREAM` is a switch of the strict set: on unless it is `0`, or `LOSAT_X_DEFAULT=0`
//! and it is not set. It acts for every TBLASTN search with `-db`.
//!
//! `LOSAT_X_TBNDBSTREAMSHADOW=1` (shadow mode, with the switch on): `read_stores` makes the copy
//! too; the check's lengths and error, and every read sequence, are compared with it. The first
//! difference stops the run (panic) after printing the line that is otherwise printed at exit:
//! `[X_SHADOW] TBNDBSTREAM sequences= reads= bases= identical`.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::OnceLock;

use anyhow::Result;

use crate::algorithm::common::database as shared;
use crate::seqdb::SeqDb;
use crate::sequence::packed_subject::PackedNucleotide;

/// The program that the reader's errors name (`database::PROGRAM`).
const PROGRAM: &str = "TBLASTN";

// No NCBI counterpart: reads the LOSAT_X_TBNDBSTREAM variable once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): a switch of the strict set (`utils/xswitch.rs`).
pub(super) fn x_tbn_db_stream() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNDBSTREAM"))
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAMSHADOW): shadow mode, with `LOSAT_X_TBNDBSTREAM` on.
pub(super) fn x_shadow() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        x_tbn_db_stream()
            && std::env::var("LOSAT_X_TBNDBSTREAMSHADOW").is_ok_and(|value| value == "1")
    })
}

static SHADOW_RAN: AtomicBool = AtomicBool::new(false);
static SHADOW_SEQUENCES: AtomicU64 = AtomicU64::new(0);
static SHADOW_READS: AtomicU64 = AtomicU64::new(0);
static SHADOW_BASES: AtomicU64 = AtomicU64::new(0);

/// The shadow line, printed at exit (`x_print_shadow_summary`) or before the panic of a difference.
fn shadow_line() -> String {
    format!(
        "[X_SHADOW] TBNDBSTREAM sequences={} reads={} bases={} identical",
        SHADOW_SEQUENCES.load(Ordering::Relaxed),
        SHADOW_READS.load(Ordering::Relaxed),
        SHADOW_BASES.load(Ordering::Relaxed),
    )
}

/// Stops the run at a difference of the shadow mode.
fn shadow_assert(same: bool, what: impl FnOnce() -> String) {
    if !same {
        eprintln!("{}", shadow_line());
        panic!("LOSAT_X_TBNDBSTREAMSHADOW: {}", what());
    }
}

/// Prints the shadow line once, if the shadow mode compared anything in this run.
pub(super) fn x_print_shadow_summary() {
    if SHADOW_RAN.swap(false, Ordering::Relaxed) {
        eprintln!("{}", shadow_line());
    }
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the database of a `-db` search, read where it is searched.
pub(super) struct XDbStream<'a> {
    db: &'a SeqDb,
    /// The length of every sequence, by ordinal id (`x_check_nucleotide_stores`).
    lengths: Vec<u32>,
    /// The copy of the shadow mode.
    copy: Option<Vec<PackedNucleotide>>,
}

impl<'a> XDbStream<'a> {
    /// Checks `db` as `read_nucleotide_stores` reads it and keeps the lengths; in the shadow mode
    /// also makes the copy (`read_stores`) and compares the check's result with it.
    pub(super) fn new(
        db: &'a SeqDb,
        read_stores: impl FnOnce() -> Result<Vec<PackedNucleotide>>,
    ) -> Result<Self> {
        let checked = shared::x_check_nucleotide_stores(db, PROGRAM);
        let copy = if x_shadow() {
            SHADOW_RAN.store(true, Ordering::Relaxed);
            let copy = read_stores();
            match (&checked, &copy) {
                (Ok(lengths), Ok(copy)) => {
                    shadow_assert(
                        lengths.len() == copy.len()
                            && lengths
                                .iter()
                                .zip(copy)
                                .all(|(&length, store)| length as usize == store.len()),
                        || "the lengths of the check differ from the copy's".to_string(),
                    );
                    SHADOW_SEQUENCES.store(lengths.len() as u64, Ordering::Relaxed);
                }
                (Err(check), Err(copy)) => shadow_assert(
                    format!("{check:#}") == format!("{copy:#}"),
                    || format!("the check's error ({check:#}) differs from the copy's ({copy:#})"),
                ),
                _ => shadow_assert(false, || {
                    "the check and the copy differ in success".to_string()
                }),
            }
            copy.ok()
        } else {
            None
        };
        Ok(Self {
            db,
            lengths: checked?,
            copy,
        })
    }

    /// The number of sequences.
    pub(super) fn count(&self) -> usize {
        self.lengths.len()
    }

    /// The number of bases of sequence `oid`.
    pub(super) fn length(&self, oid: usize) -> usize {
        self.lengths[oid] as usize
    }

    /// Sequence `oid`, packed as the copy packs it, with one read.
    pub(super) fn read(&self, oid: usize) -> Result<PackedNucleotide> {
        let ordinal = u32::try_from(oid).expect("a database ordinal id fits Uint4");
        let store = shared::x_read_nucleotide_store(self.db, ordinal, PROGRAM)?;
        if let Some(copy) = &self.copy {
            shadow_assert(store == copy[oid], || {
                format!("database sequence {oid} differs from the copy")
            });
            SHADOW_READS.fetch_add(1, Ordering::Relaxed);
            SHADOW_BASES.fetch_add(store.len() as u64, Ordering::Relaxed);
        }
        Ok(store)
    }
}
