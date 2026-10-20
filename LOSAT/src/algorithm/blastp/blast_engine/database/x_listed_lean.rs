//! EXPERIMENT (LOSAT_X_BPLISTEDLEAN / LOSAT_X_BPLISTEDLEANSHADOW): the reports of a `-db` search
//! read the letters of the listed subjects from their stored residues and hold the titles of the
//! reported subjects only.
//!
//! Without the switch, `ListedSubjects::read_titles` makes one `FastaRecord` per listed subject
//! (every subject of the preliminary hit lists): its letters in ASCII and its stored title, empty
//! for a subject in no final hit list. The reports then copy every title into an `Arc` twice, once
//! for outfmt 0 (`subject_titles`) and once for the `stitle` field (`ListedTitles::stitles`). With
//! the UniRef50 database, the first query batch (91 queries) of 100 at 32 threads lists 70,701
//! subjects of 50.8 MB of letters: the records added 41 MiB and the two title copies 10 MiB of
//! anonymous memory on top of the listed sequences (`ListedSubjects::encoded`, which the redo reads and the reports keep), all in
//! the main heap in about 0.2 s (v4-D12 s3 probe, workstation).
//!
//! With the switch on, `ListedSubjects::x_read_titles_lean` reads the same titles, in the same
//! order and with the same errors (the subjects of the final hit lists, query by query, each once),
//! into one `Arc` per reported subject; the other subjects share one empty title. The reports
//! borrow these titles for outfmt 0 and for `stitle`, and read the letters of a subject from its
//! stored residues by `ncbistdaa_to_ascii` (`blast_engine.rs` `ShownSubjects`), which is how
//! `read_titles` makes the records' letters. Only two readers take letters from the records: the
//! recount of the HSPs with pyrrolysine (`recount_pyrrolysine_identities`) and the O of the rendered
//! rows (`restore_displayed_o`), both for the subjects of the final hit lists. Every other reader
//! of a listed subject takes the stored residues (`subjects`) or its title by rank. So every letter
//! and title that a report reads is the byte that the records held, and the reports do not change.
//!
//! No NCBI counterpart for the holding: NCBI reads each sequence and definition line from the
//! mapped database files when it formats a hit (`CSeqDBVol::x_GetHdrAsn1Binary`, quoted above
//! `read_titles`).
//!
//! `LOSAT_X_BPLISTEDLEAN` is a switch of the strict set: on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set. It acts for every task of BLASTP with `-db`, in the
//! search of all queries and in the blocks of `LOSAT_X_BPQBLOCK` (not with ABI v1's titles, which
//! read the records).
//!
//! `LOSAT_X_BPLISTEDLEANSHADOW=1` (shadow mode, implies the switch): `read_titles` runs too, first,
//! and its records are compared with the lean titles rank by rank: every title, every flag of
//! replaced title bytes, and every letter against the stored residue's `ncbistdaa_to_ascii`. The
//! first difference stops the run (panic) after printing the line that is otherwise printed at
//! exit: `[X_SHADOW] BPLISTEDLEAN searches= ranks= reported= letters= identical` (`ranks`: listed
//! subjects compared; `reported`: those with a title read; `letters`: letters compared).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};

use anyhow::Result;

use super::{
    ncbistdaa_to_ascii, shared, BlastpDatabase, BlastpHitList, ListedShown, ListedSubjects,
    ListedTitles, PROGRAM,
};

/// EXPERIMENT (LOSAT_X_BPLISTEDLEAN): how the reports of a `-db` search hold the listed subjects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::algorithm::blastp::blast_engine) enum XListedLean {
    /// One record per listed subject, with its letters (`ListedSubjects::read_titles`).
    Off,
    /// The titles of the reported subjects; the letters from the stored residues.
    On,
    /// `On`, with the records of `Off` made too and compared.
    Shadow,
}

// No NCBI counterpart: reads the LOSAT_X_BPLISTEDLEAN / LOSAT_X_BPLISTEDLEANSHADOW switches once;
// it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLISTEDLEAN): a strict-set switch, on unless it is `0`, or
/// `LOSAT_X_DEFAULT=0` and it is not set; `LOSAT_X_BPLISTEDLEANSHADOW` (any value) selects the
/// shadow mode.
pub(in crate::algorithm::blastp::blast_engine) fn x_bp_listed_lean() -> XListedLean {
    static MODE: OnceLock<XListedLean> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BPLISTEDLEANSHADOW").is_some() {
            XListedLean::Shadow
        } else if crate::utils::xswitch::on("LOSAT_X_BPLISTEDLEAN") {
            XListedLean::On
        } else {
            XListedLean::Off
        }
    })
}

impl ListedSubjects {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1582 and
    // c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2302-2310, quoted above `read_titles`:
    // NCBI reads the definition lines of the final hit lists, query by query, in the order of each
    // hit list; these are the same reads.
    /// EXPERIMENT (LOSAT_X_BPLISTEDLEAN): `read_titles` without the records: the definition lines
    /// of the subjects of the final hit lists (`hit_lists`), read in the same order with the same
    /// errors, one `Arc` each; the other subjects share one empty title.
    pub(in crate::algorithm::blastp::blast_engine) fn x_read_titles_lean(
        &self,
        database: &BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
    ) -> Result<ListedTitles<'_>> {
        let ranks = self.ordinals.len();
        let mut titles: Vec<Arc<[u8]>> = vec![Arc::from(&b""[..]); ranks];
        let mut replaced_title_bytes = vec![false; ranks];
        let mut read = vec![false; ranks];
        for (query, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            for hsp_list in &hit_list.hsplist_array[..count] {
                let rank = hsp_list.oid as usize;
                if std::mem::replace(&mut read[rank], true) {
                    continue;
                }
                let stored = shared::read_title(&database.db, self.ordinals[rank], query, PROGRAM)?;
                titles[rank] = Arc::from(stored.title);
                replaced_title_bytes[rank] = stored.replaced_bytes;
            }
        }
        if x_bp_listed_lean() == XListedLean::Shadow {
            X_REPORTED.fetch_add(
                read.iter().filter(|&&read| read).count() as u64,
                Ordering::Relaxed,
            );
        }
        Ok(ListedTitles {
            listed: self,
            shown: ListedShown::Lean(titles),
            replaced_title_bytes,
        })
    }
}

static X_SEARCHES: AtomicU64 = AtomicU64::new(0);
static X_RANKS: AtomicU64 = AtomicU64::new(0);
static X_REPORTED: AtomicU64 = AtomicU64::new(0);
static X_LETTERS: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: prints counters and changes no value NCBI computes.
fn x_shadow_print_line(identical: bool) {
    eprintln!(
        "[X_SHADOW] BPLISTEDLEAN searches={} ranks={} reported={} letters={}{}",
        X_SEARCHES.load(Ordering::Relaxed),
        X_RANKS.load(Ordering::Relaxed),
        X_REPORTED.load(Ordering::Relaxed),
        X_LETTERS.load(Ordering::Relaxed),
        if identical { " identical" } else { "" }
    );
}

// No NCBI counterpart: prints counters and changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLISTEDLEANSHADOW): the shadow mode's line at exit.
pub(in crate::algorithm::blastp) fn x_bp_listed_lean_print_summary() {
    if x_bp_listed_lean() == XListedLean::Shadow {
        x_shadow_print_line(true);
    }
}

// No NCBI counterpart: the shadow mode's check stops the run on a difference.
fn x_shadow_difference(message: String) -> ! {
    // The panic ends the run before `main` prints the line at exit: it is printed here.
    x_shadow_print_line(false);
    panic!("[X_SHADOW] BPLISTEDLEAN {message}");
}

// No NCBI counterpart: a check of the shadow mode; it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLISTEDLEANSHADOW): the lean titles (`lean`) against the records of
/// `read_titles` (`records`) of the same search, rank by rank: the title, the flag of replaced
/// title bytes and every letter (the stored residue's `ncbistdaa_to_ascii`).
pub(in crate::algorithm::blastp::blast_engine) fn x_shadow_compare(
    lean: &ListedTitles<'_>,
    records: &ListedTitles<'_>,
) {
    X_SEARCHES.fetch_add(1, Ordering::Relaxed);
    let listed = lean.listed;
    let (ListedShown::Lean(titles), ListedShown::Records(records_of)) =
        (&lean.shown, &records.shown)
    else {
        x_shadow_difference("the lean titles against the records: wrong kinds".to_string());
    };
    if titles.len() != listed.ordinals.len() || records_of.len() != listed.ordinals.len() {
        x_shadow_difference(format!(
            "ranks: lean titles {} / records {} / listed {}",
            titles.len(),
            records_of.len(),
            listed.ordinals.len()
        ));
    }
    for (rank, record) in records_of.iter().enumerate() {
        let ordinal = listed.ordinals[rank];
        if titles[rank][..] != record.title[..] {
            x_shadow_difference(format!("rank {rank} (ordinal {ordinal}): title"));
        }
        if lean.replaced_title_bytes[rank] != records.replaced_title_bytes[rank] {
            x_shadow_difference(format!(
                "rank {rank} (ordinal {ordinal}): replaced title bytes"
            ));
        }
        let stored = &listed.encoded[rank];
        let residues = &stored.aa_seq[1..=stored.aa_len];
        let letters = record.seq();
        if residues.len() != letters.len()
            || residues
                .iter()
                .zip(letters)
                .any(|(&residue, &letter)| ncbistdaa_to_ascii(residue) as u8 != letter)
        {
            x_shadow_difference(format!("rank {rank} (ordinal {ordinal}): letters"));
        }
        X_LETTERS.fetch_add(letters.len() as u64, Ordering::Relaxed);
    }
    X_RANKS.fetch_add(records_of.len() as u64, Ordering::Relaxed);
}
