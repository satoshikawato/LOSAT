//! EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the reports of a `-db` search read the letters of the
//! listed subjects from their stored residues and hold the titles of the reported subjects only.
//!
//! Without the switch, `ListedSubjects::read_titles` makes one `FastaRecord` per listed subject
//! (every subject of the preliminary hit lists): its letters in ASCII and its stored title, empty
//! for a subject in no final hit list. The reports then copy every title into an `Arc` twice, once
//! for outfmt 0 (`subject_titles`) and once for the `stitle` field (`ListedTitles::stitles`). With
//! the UniRef50 database, the first query batch (91 queries) of 100 at 32 threads lists 70,701
//! subjects of 50.8 MB of letters: the records added 41 MiB and the two title copies 10 MiB of
//! anonymous memory on top of the listed sequences (`ListedSubjects::encoded`, which the redo read
//! and the reports kept before `LOSAT_X_BPREDODBSUBJ` graduated), all in the main heap in about 0.2 s (v4-D12
//! s3 probe, workstation).
//!
//! With the switch on, `ListedSubjects::x_read_titles_lean` reads the same titles, in the same
//! order and with the same errors (the subjects of the final hit lists, query by query, each once),
//! into one `Arc` per reported subject; the other subjects share one empty title. The reports
//! borrow these titles for outfmt 0 and for `stitle`, and read the letters of a subject from its
//! stored residues by `ncbistdaa_to_ascii` (`blast_engine.rs` `ShownSubjects`), which is how
//! `read_titles` makes the records' letters. Only two readers take letters from the records: the
//! recount of the HSPs with pyrrolysine or a gap character (`recount_shown_identities`) and the O
//! of the rendered rows (`restore_displayed_o`), both for the subjects of the final hit lists.
//! Every other reader of a listed subject takes the stored residues (`subjects`) or its title by
//! rank. So every letter and title that a report reads is the byte that the records held, and the
//! reports do not change.
//!
//! No NCBI counterpart for the holding: NCBI reads each sequence and definition line from the
//! mapped database files when it formats a hit (`CSeqDBVol::x_GetHdrAsn1Binary`, quoted above
//! `read_titles`).
//!
//! `LOSAT_X_BPLISTEDLEAN` is a switch of the strict set: on unless it is `0`.
//! It acts for every task of BLASTP with `-db`, in the
//! search of all queries and in the blocks of `LOSAT_X_BPQBLOCK` (not with ABI v1's titles, which
//! read the records).

use std::sync::{Arc, OnceLock};

use anyhow::Result;

#[cfg(test)]
use super::ncbistdaa_to_ascii;
use super::{
    BlastpDatabase, BlastpHitList, ListedShown, ListedSubjects, ListedTitles, ParsedSubjects,
};

// No NCBI counterpart: reads the LOSAT_X_BPLISTEDLEAN switch once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLISTEDLEAN): a strict-set switch, on unless it is `0`.
pub(in crate::algorithm::blastp::blast_engine) fn x_bp_listed_lean() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPLISTEDLEAN"))
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
        let mut parsed = ParsedSubjects::default();
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
                let (title, replaced) =
                    parsed.read(database, rank, self.ordinals[rank], query, ranks)?;
                titles[rank] = Arc::from(title);
                replaced_title_bytes[rank] = replaced;
            }
        }
        Ok(ListedTitles {
            listed: self,
            shown: ListedShown::Lean(titles),
            replaced_title_bytes,
            parsed: parsed.by_rank,
        })
    }
}

// No NCBI counterpart: a check of the unit tests; it changes no value NCBI computes.
/// The lean titles (`lean`) against the records of `read_titles` (`records`) of the same search,
/// rank by rank: the title, the flag of replaced title bytes and every letter (the stored
/// residue's `ncbistdaa_to_ascii`).
#[cfg(test)]
pub(in crate::algorithm::blastp::blast_engine) fn x_assert_same_titles_and_letters(
    lean: &ListedTitles<'_>,
    records: &ListedTitles<'_>,
) {
    let listed = lean.listed;
    let (ListedShown::Lean(titles), ListedShown::Records(records_of)) =
        (&lean.shown, &records.shown)
    else {
        panic!("the lean titles against the records: wrong kinds");
    };
    assert!(
        titles.len() == listed.ordinals.len() && records_of.len() == listed.ordinals.len(),
        "ranks: lean titles {} / records {} / listed {}",
        titles.len(),
        records_of.len(),
        listed.ordinals.len()
    );
    for (rank, record) in records_of.iter().enumerate() {
        let ordinal = listed.ordinals[rank];
        assert!(
            titles[rank][..] == record.title[..],
            "rank {rank} (ordinal {ordinal}): title"
        );
        assert_eq!(
            lean.replaced_title_bytes[rank], records.replaced_title_bytes[rank],
            "rank {rank} (ordinal {ordinal}): replaced title bytes"
        );
        let stored = &listed.encoded[rank];
        let residues = &stored.aa_seq[1..=stored.aa_len];
        let letters = record.seq();
        assert!(
            residues.len() == letters.len()
                && residues
                    .iter()
                    .zip(letters)
                    .all(|(&residue, &letter)| ncbistdaa_to_ascii(residue) as u8 == letter),
            "rank {rank} (ordinal {ordinal}): letters"
        );
    }
}
