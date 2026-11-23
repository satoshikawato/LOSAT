//! The listed subjects of a `-db` search without a copy of their residues; the redo reads each
//! match's subject from the database (`blastp/x_redo_subjects.rs`).

use std::sync::atomic::Ordering;

use anyhow::Result;

use super::{
    encoded_stored, listed_ordinals, rename_to_ranks, BlastpDatabase, BlastpHitList,
    EncodedProtein, ListedSubjects, XPrelimHsp,
};
use crate::algorithm::blastp::x_redo_subjects::{XDbSubjects, XRedoSubjects};
use crate::algorithm::common::database::database_error;

/// The residues that the reports of a search read
/// (`ListedSubjects::x_read_shown`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::algorithm::blastp::blast_engine) enum XShownResidues {
    /// Those of the subjects of the final hit lists with a pyrrolysine (O) or a gap character, and
    /// of every subject of an HSP whose query has an O: the reports read the letters only to
    /// recount the HSPs whose query has an O or whose subject has an O or a gap character
    /// (`recount_shown_identities`; the tabular rows of the default fields, the Seq-aligns).
    WithO,
    /// Those of every subject of the final hit lists: a report renders the alignments.
    Final,
    /// Those of every listed subject: the reports make records (`ListedSubjects::read_titles`).
    Listed,
}

impl ListedSubjects {
    /// The ordinal ids of the subjects of `hit_lists` and `x_compact_hit_lists` (the lists not
    /// expanded yet, `x_prelim_compact.rs`) and their lengths, read from the index, without a copy
    /// of the residues; the subjects of both are renamed to their ranks. The subjects that NCBI's
    /// fetches of the redo or the traceback exclude leave both lists first, each list leaving a
    /// mark of its place (`exclude_fetched`, `fetch_excluded`; a header that the fetch cannot read
    /// is LOSAT's refusal where the redo or the traceback fetches it, `note_fetched`).
    pub(in crate::algorithm::blastp::blast_engine) fn x_list(
        database: &BlastpDatabase,
        hit_lists: &mut [Option<BlastpHitList>],
        x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    ) -> Result<Self> {
        let (fetch_excluded, fetched_utf8) =
            super::exclude_fetched(database, hit_lists, x_compact_hit_lists);
        let ordinals = listed_ordinals(hit_lists, x_compact_hit_lists);
        let x_lengths = ordinals
            .iter()
            .map(|&ordinal| {
                database
                    .db
                    .protein_length(ordinal)
                    .map_err(|error| database_error(error, database.program))
            })
            .collect::<Result<Vec<u32>>>()?;
        rename_to_ranks(&ordinals, hit_lists, x_compact_hit_lists);
        Ok(Self {
            x_recounted: ordinals.iter().map(|_| Default::default()).collect(),
            ordinals,
            encoded: Vec::new(),
            x_lengths,
            fetched_utf8,
            earlier_utf8: Default::default(),
            fetch_excluded,
        })
    }

    /// The subjects as the redo reads them: the database by rank.
    /// `hard_mask`: the algorithm of `-db_hard_mask`, whose masks the redo's reads write `X`
    /// over (`XDbSubjects::mask_hard`).
    pub(in crate::algorithm::blastp::blast_engine) fn redo_subjects<'a>(
        &'a self,
        database: &'a BlastpDatabase,
        hard_mask: Option<i32>,
    ) -> XRedoSubjects<'a> {
        XRedoSubjects::Database(XDbSubjects {
            db: &database.db,
            ordinals: &self.ordinals,
            lengths: &self.x_lengths,
            recounted: &self.x_recounted,
            program: database.program,
            hard_mask,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2302-2310, quoted
    // above `read_titles`: NCBI reads the database again when it formats the results of a batch.
    /// After the redo, `encoded` for the reports: every listed subject with its length, and with
    /// its residues read from the database where the reports read them (`residues`;
    /// `query_recounted`: whether the recount counts each query's HSPs again, by query index,
    /// `ShownLetters::recounted`, read with `WithO`).
    pub(in crate::algorithm::blastp::blast_engine) fn x_read_shown(
        &mut self,
        database: &BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
        residues: XShownResidues,
        query_recounted: &[bool],
    ) -> Result<()> {
        let ranks = self.ordinals.len();
        let mut reported = vec![residues == XShownResidues::Listed; ranks];
        // With `WithO`, the subjects of the HSPs of a query that the recount counts again (an O, a
        // gap residue, Bioseq letters): it reads all their letters.
        let mut kept = vec![residues != XShownResidues::WithO; ranks];
        for hit_list in hit_lists.iter().flatten() {
            for hsp_list in &hit_list.hsplist_array {
                reported[hsp_list.oid as usize] = true;
                for hsp in &hsp_list.hsps {
                    reported[hsp.s_idx as usize] = true;
                    if query_recounted.get(hsp.q_idx as usize) == Some(&true) {
                        kept[hsp.s_idx as usize] = true;
                    }
                }
            }
        }
        let mut encoded = Vec::with_capacity(ranks);
        for (rank, (&ordinal, &length)) in self.ordinals.iter().zip(&self.x_lengths).enumerate() {
            let length_only = EncodedProtein {
                aa_seq: Vec::new(),
                aa_len: length as usize,
            };
            if !reported[rank] {
                encoded.push(length_only);
                continue;
            }
            // With `WithO`, a subject that no read of the redo found an O or a gap character in (the
            // redo read every subject of the final hit lists) and no HSP of a recounted query names
            // keeps its length only, without a read.
            if !kept[rank] && !self.x_recounted[rank].load(Ordering::Relaxed) {
                encoded.push(length_only);
                continue;
            }
            let stored = database
                .db
                .read_sequence(ordinal)
                .map_err(|error| database_error(error, database.program))?;
            let stored = stored.protein(0);
            encoded.push(encoded_stored(stored));
        }
        self.encoded = encoded;
        Ok(())
    }
}
