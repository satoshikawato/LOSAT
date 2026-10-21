//! EXPERIMENT (LOSAT_X_BPREDODBSUBJ / LOSAT_X_BPREDODBSUBJSHADOW): the listed subjects of a `-db`
//! search without the copy of their residues; the redo reads each match's subject from the
//! database (`blastp/x_redo_subjects.rs`, which describes the switch).

use anyhow::Result;

use super::{
    database_error, encoded_stored, listed_ordinals, read_listed, rename_to_ranks, BlastpDatabase,
    BlastpHitList, EncodedProtein, ListedSubjects, XPrelimHsp, PROGRAM,
};
use crate::algorithm::blastp::x_redo_subjects::{
    x_has_stored_o, x_shadow_compare_lengths, x_shadow_compare_residues, XDbSubjects, XRedoDbSubj,
    XRedoSubjects,
};

/// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): the residues that the reports of a search read
/// (`ListedSubjects::x_read_shown`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::algorithm::blastp::blast_engine) enum XShownResidues {
    /// Those of the subjects of the final hit lists with a pyrrolysine (O), and of every subject of
    /// an HSP whose query has an O: the reports read the letters only to recount the HSPs whose
    /// query or subject has an O (`recount_pyrrolysine_identities`; the tabular rows of the
    /// default fields, the Seq-aligns).
    WithO,
    /// Those of every subject of the final hit lists: a report renders the alignments.
    Final,
    /// Those of every listed subject: the reports make records (`ListedSubjects::read_titles`).
    Listed,
}

impl ListedSubjects {
    /// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): `read` without the copy of the residues: the ordinal ids
    /// of the subjects of `hit_lists` and `x_compact_hit_lists` and their lengths, read from the
    /// index; the subjects of both are renamed to their ranks. The shadow mode (`mode`) makes the
    /// copy of `read` too, first, and compares the lengths with it.
    pub(in crate::algorithm::blastp::blast_engine) fn x_list(
        database: &BlastpDatabase,
        hit_lists: &mut [Option<BlastpHitList>],
        x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
        mode: XRedoDbSubj,
    ) -> Result<Self> {
        let ordinals = listed_ordinals(hit_lists, x_compact_hit_lists);
        let x_shadow_copy = if mode == XRedoDbSubj::Shadow {
            read_listed(database, &ordinals)?
        } else {
            Vec::new()
        };
        let x_lengths = ordinals
            .iter()
            .map(|&ordinal| database.db.protein_length(ordinal).map_err(database_error))
            .collect::<Result<Vec<u32>>>()?;
        if mode == XRedoDbSubj::Shadow {
            x_shadow_compare_lengths(&ordinals, &x_lengths, &x_shadow_copy);
        }
        rename_to_ranks(&ordinals, hit_lists, x_compact_hit_lists);
        Ok(Self {
            ordinals,
            encoded: Vec::new(),
            x_mode: mode,
            x_lengths,
            x_shadow_copy,
        })
    }

    /// The subjects as the redo reads them: the copy (`read`), or, with LOSAT_X_BPREDODBSUBJ
    /// (`x_list`), the database by rank.
    pub(in crate::algorithm::blastp::blast_engine) fn redo_subjects<'a>(
        &'a self,
        database: &'a BlastpDatabase,
    ) -> XRedoSubjects<'a> {
        match self.x_mode {
            XRedoDbSubj::Off => XRedoSubjects::Encoded(&self.encoded),
            XRedoDbSubj::On | XRedoDbSubj::Shadow => XRedoSubjects::Database(XDbSubjects {
                db: &database.db,
                ordinals: &self.ordinals,
                lengths: &self.x_lengths,
                program: PROGRAM,
                shadow: (self.x_mode == XRedoDbSubj::Shadow).then_some(&self.x_shadow_copy[..]),
            }),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2302-2310, quoted
    // above `read_titles`: NCBI reads the database again when it formats the results of a batch.
    /// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): after the redo, `encoded` for the reports (nothing to do
    /// without the switch): every listed subject with its length, and with its residues read from
    /// the database where the reports read them (`residues`; `query_has_o`: whether each query's
    /// record has an O, by query index, read with `WithO`). The shadow mode compares each read with
    /// the copy.
    pub(in crate::algorithm::blastp::blast_engine) fn x_read_shown(
        &mut self,
        database: &BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
        residues: XShownResidues,
        query_has_o: &[bool],
    ) -> Result<()> {
        if self.x_mode == XRedoDbSubj::Off {
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
        let shadow = self.x_mode == XRedoDbSubj::Shadow;
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
            let stored = database.db.read_sequence(ordinal).map_err(database_error)?;
            let stored = stored.protein(0);
            if shadow {
                x_shadow_compare_residues(rank, ordinal, &self.x_shadow_copy[rank], stored);
            }
            if !kept[rank] && !x_has_stored_o(stored) {
                encoded.push(length_only);
            } else {
                encoded.push(encoded_stored(stored));
            }
        }
        self.encoded = encoded;
        Ok(())
    }
}
