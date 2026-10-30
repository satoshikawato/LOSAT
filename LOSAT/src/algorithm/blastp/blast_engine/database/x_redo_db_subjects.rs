//! EXPERIMENT (LOSAT_X_BPREDODBSUBJ): the listed subjects of a `-db` search without the copy of
//! their residues; the redo reads each match's subject from the database
//! (`blastp/x_redo_subjects.rs`, which describes the switch).

use std::sync::atomic::Ordering;

use anyhow::Result;

use super::{
    database_error, encoded_stored, listed_ordinals, rename_to_ranks, BlastpDatabase,
    BlastpHitList, EncodedProtein, ListedSubjects, XPrelimHsp, PROGRAM,
};
use crate::algorithm::blastp::x_redo_subjects::{XDbSubjects, XRedoSubjects};

/// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): the residues that the reports of a search read
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
    /// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): `read` without the copy of the residues: the ordinal ids
    /// of the subjects of `hit_lists` and `x_compact_hit_lists` and their lengths, read from the
    /// index; the subjects of both are renamed to their ranks.
    pub(in crate::algorithm::blastp::blast_engine) fn x_list(
        database: &BlastpDatabase,
        hit_lists: &mut [Option<BlastpHitList>],
        x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    ) -> Result<Self> {
        let ordinals = listed_ordinals(hit_lists, x_compact_hit_lists);
        let x_lengths = ordinals
            .iter()
            .map(|&ordinal| database.db.protein_length(ordinal).map_err(database_error))
            .collect::<Result<Vec<u32>>>()?;
        rename_to_ranks(&ordinals, hit_lists, x_compact_hit_lists);
        Ok(Self {
            x_recounted: ordinals.iter().map(|_| Default::default()).collect(),
            ordinals,
            encoded: Vec::new(),
            x_database: true,
            x_lengths,
        })
    }

    /// The subjects as the redo reads them: the copy (`read`), or, with LOSAT_X_BPREDODBSUBJ
    /// (`x_list`), the database by rank.
    pub(in crate::algorithm::blastp::blast_engine) fn redo_subjects<'a>(
        &'a self,
        database: &'a BlastpDatabase,
    ) -> XRedoSubjects<'a> {
        if self.x_database {
            XRedoSubjects::Database(XDbSubjects {
                db: &database.db,
                ordinals: &self.ordinals,
                lengths: &self.x_lengths,
                recounted: &self.x_recounted,
                program: PROGRAM,
            })
        } else {
            XRedoSubjects::Encoded(&self.encoded)
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2302-2310, quoted
    // above `read_titles`: NCBI reads the database again when it formats the results of a batch.
    /// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): after the redo, `encoded` for the reports (nothing to do
    /// without the switch): every listed subject with its length, and with its residues read from
    /// the database where the reports read them (`residues`; `query_has_o`: whether each query's
    /// record has an O, by query index, read with `WithO`).
    pub(in crate::algorithm::blastp::blast_engine) fn x_read_shown(
        &mut self,
        database: &BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
        residues: XShownResidues,
        query_has_o: &[bool],
    ) -> Result<()> {
        if !self.x_database {
            return Ok(());
        }
        let ranks = self.ordinals.len();
        let mut reported = vec![residues == XShownResidues::Listed; ranks];
        // With `WithO`, the subjects of the HSPs whose query has an O: the recount reads all their
        // letters.
        let mut kept = vec![residues != XShownResidues::WithO; ranks];
        for hit_list in hit_lists.iter().flatten() {
            for hsp_list in &hit_list.hsplist_array {
                reported[hsp_list.oid as usize] = true;
                for hsp in &hsp_list.hsps {
                    reported[hsp.s_idx as usize] = true;
                    if query_has_o.get(hsp.q_idx as usize) == Some(&true) {
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
            // redo read every subject of the final hit lists) and no HSP of a query with an O names
            // keeps its length only, without a read.
            if !kept[rank] && !self.x_recounted[rank].load(Ordering::Relaxed) {
                encoded.push(length_only);
                continue;
            }
            let stored = database.db.read_sequence(ordinal).map_err(database_error)?;
            let stored = stored.protein(0);
            encoded.push(encoded_stored(stored));
        }
        self.encoded = encoded;
        Ok(())
    }
}
