//! The volumes of a profile database as `rpsblast` and `rpstblastn` search them (series PR,
//! items PS11b and PS12a): NCBI's order of the volumes and their first ordinal ids, the search
//! units and their statistics under the alias files' keys (`RpsUnits`, line RG D1a-4), the set-up
//! of one volume (its files, LOSAT's refusals of what NCBI's search does not reproducibly run, the
//! `.psq` letters of its profiles), the masking of the letters that a profile has no column
//! for, the merge of the volumes' results and their ordinal ids for the reports. One owner for
//! both programs (review A of PS12a, N3).

use anyhow::{bail, Result};

use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::seqdb::rps::{canonical_volume_path, RpsOpenError, RpsVolumeFault, RpsVolumeFiles};

/// NCBI's own errors of a volume's files, which its search of a query batch raises (`BLAST
/// engine error`, exit status 3).
pub(crate) enum RpsVolumeError {
    /// NCBI's bytes of a missing, short or unreadable file (its path as the system's bytes).
    Files(Vec<u8>),
    /// An `.aux` without a matrix name (`RpsVolumeFault::NoMatrix`): NCBI's score block fails for
    /// each query of the batch ("Unknown error code -1").
    NoMatrix,
}
use crate::seqdb::{SeqDb, Volume};

/// The volumes of a profile database in the order of their search.
pub(crate) struct RpsVolumeSet<'a> {
    /// The database's volumes in its own order (that of the ordinal ids).
    pub(crate) volumes: &'a [Volume],
    /// The first ordinal id of each volume.
    pub(crate) starts: Vec<u32>,
    /// The volumes in the order NCBI searches them.
    pub(crate) order: Vec<usize>,
}

impl<'a> RpsVolumeSet<'a> {
    pub(crate) fn of(db: &'a SeqDb) -> Self {
        let db_volumes = db.volumes();
        // The volumes in the order of their paths as strings (a name given twice is searched
        // once).
        let mut order: Vec<usize> = (0..db_volumes.len()).collect();
        order.sort_by(|&a, &b| {
            crate::seqdb::SeqDb::compare_volume_paths(db_volumes[a].base(), db_volumes[b].base())
        });
        Self::with_order(db, order)
    }

    /// The set of `db`'s volumes in the search order `order` that `of` gave them before
    /// (`RpsUnits::order`), without sorting them again.
    pub(crate) fn with_order(db: &'a SeqDb, order: Vec<usize>) -> Self {
        let db_volumes = db.volumes();
        // The first ordinal id of each volume (the database's order).
        let mut starts = Vec::with_capacity(db_volumes.len());
        let mut next_start = 0u32;
        for volume in db_volumes {
            starts.push(next_start);
            next_start += volume.num_oids();
        }
        Self {
            volumes: db_volumes,
            starts,
            order,
        }
    }
    /// The first volume searched (in path order), whose statistics the merged result keeps and
    /// the reports print with `-dbsize` on several volumes (`s_CombineSearchSets`'s
    /// `thread_results[0]`).
    pub(crate) fn first(&self) -> &'a Volume {
        &self.volumes[self.order[0]]
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:905-937 (cited at
// `SeqDb::statistics_total_length`) and c++/src/algo/blast/core/blast_engine.c:1129-1143
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:131-139 (a bare
// volume walks no alias file: its index files store no shortest length, and it has no STATS key)
/// The numbers of one search unit's `BlastSeqSrc` that NCBI's statistics and RPS engine read
/// (`RpsUnits`; line RG D1a-4, model M2), with NCBI's types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RpsUnit {
    /// `BlastSeqSrcGetTotLenStats`, else `BlastSeqSrcGetTotLen`, as NCBI's `Int8` before the
    /// fallback to OID 0 (`T0`): the Gumbel block's length and the search strategy's.
    pub(crate) total: i64,
    /// `BLAST_GapAlignSetUp`'s total and number of sequences (`T`, `n`): `total` and
    /// `STATS_NSEQ`, else `GetNumSeqs`, for a positive `total`, else the length of OID 0 and 1
    /// (0 and 0 without an OID).
    pub(crate) setup_total: i64,
    pub(crate) setup_count: i32,
    /// `STATS_NSEQ`, else `GetNumSeqs`, without the fallback (`SplitQuery_SetEffectiveSearchSpace`).
    pub(crate) request_count: i32,
    /// `GetMinLength`: the subject length of the hit cutoff with the Gumbel block.
    pub(crate) min_length: i32,
    /// The RPS engine's `dbsize`, `BlastSeqSrcGetTotLen` plus `BlastSeqSrcGetNumSeqs`, widened from
    /// NCBI's `Int8` (rule P2: NCBI's sum overflows for an alias file's `LENGTH` of 2^63 - 1 and
    /// an `NSEQ` of 1 or more, undefined behaviour; widened it is above `INT4_MAX`, NCBI's -3 and
    /// LOSAT's refusal of the double free, `RpsUnit::check_engine_length`), whose `Int4` is the
    /// concatenated length `Lc` (`concat_length`).
    pub(crate) engine_size: i128,
    /// Whether the unit's `BlastSeqSrcGetTotLen` (`GetTotalLength`, no STATS key) is above 0, the
    /// gate of the sort after the traceback's filters (`RpsRun::sorts_results`): not for an alias
    /// file's `LENGTH` of 0 or of 2^63 or more over one volume.
    pub(crate) sorts_results: bool,
}

impl RpsUnit {
    /// The concatenated length `Lc` of the unit's search, `concat_db.length` set from `dbsize` cast
    /// to `Int4` (blast_engine.c:1143, cited above): the low 32 bits of `engine_size`. A bare
    /// volume's is the letters of its index and one delimiter row per profile, the rows of its
    /// `.rps` file; an alias file's `LENGTH` and `NSEQ` over one volume move it below the rows or
    /// above them, also to a negative value. The RPS engine reads it in three places
    /// (`RpsRun::concat_length`): the size of the diagonal table, the bound of the ungapped
    /// extensions to the right and the root of the gapped stage's interval tree.
    pub(crate) fn concat_length(&self) -> i32 {
        self.engine_size as i32
    }

    /// The unit of the `-db` name with all its alias keys, the only unit of one volume
    /// (`SeqDb::statistics`).
    fn of_name(db: &SeqDb) -> Self {
        let statistics = db.statistics();
        Self {
            total: statistics.gumbel_db_length,
            setup_total: statistics.setup_total,
            setup_count: statistics.setup_count,
            request_count: statistics.request_count,
            min_length: statistics.min_length,
            engine_size: i128::from(statistics.engine_db_length) + i128::from(db.ncbi_num_seqs()),
            sorts_results: statistics.sorts_results(),
        }
    }

    /// The unit of a bare volume of several (`s_ModifyVolumePaths`): its own totals, no STATS
    /// key, the shortest length `BLAST_SEQSRC_MINLENGTH`.
    fn of_volume(volume: &Volume) -> std::result::Result<Self, crate::seqdb::SeqDbError> {
        let total = volume.index().volume_length as i64;
        let count = volume.num_oids() as i32;
        let (setup_total, setup_count) = if total > 0 {
            (total, count)
        } else if volume.num_oids() > 0 {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:932-938
            let len = volume.protein_length_i32(0)?;
            (if len < 0 { -1 } else { i64::from(len) }, 1)
        } else {
            (0, 0)
        };
        Ok(Self {
            total,
            setup_total,
            setup_count,
            request_count: count,
            min_length: crate::seqdb::BLAST_SEQSRC_MINLENGTH,
            engine_size: i128::from(total) + i128::from(count),
            sorts_results: total > 0,
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:250-262, 380-399
// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:349-371 (cited above
// `search_sets`): `x_AdjustDbSize`, before the volumes of several.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-727 (cited at
// `stats::setup::eff_lengths_database`) and 912-923 (the Gumbel block's `db_length`, cited at
// `blastx::parameters::gumbel_db_length`)
/// The search units of an RPS search (`rpsblast`, `rpstblastn`; line RG D1a-4, model M1) and the
/// statistics of each, the one owner of both programs: NCBI's `CLocalRPSBlast` searches the
/// `-db` name itself when it has one volume (every alias key reaches the search), else each of
/// the volumes bare, in path order (`RpsVolumeSet::order`: their own totals, no STATS key, the
/// shortest length 10, no OID range), with the database-wide length `D` and count `M` of
/// `x_AdjustDbSize` in the options (several volumes, neither `-searchsp` nor `-dbsize`). Each
/// unit's search reads its own numbers under those options: the effective lengths
/// (`eff_lengths`), the Gumbel block's length (`gumbel_length`), the shortest length; the first
/// unit's statistics are those of the merged result (the footers, `thread_results[0]`), the
/// last unit's search space is what the options handle keeps for the search strategy
/// (`request_space`). Built once per run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RpsUnits {
    /// The units in the order of their search.
    units: Vec<RpsUnit>,
    /// `x_AdjustDbSize`'s `D` (`STATS_TOTLEN`, else `GetTotalLength`, NCBI's `Uint8` as the
    /// option's `Int8`) and `M` (`STATS_NSEQ`, else `GetNumSeqs`), where it runs.
    adjusted: Option<(i64, i32)>,
    /// The `DbLength` option of the units' searches: `D`, else `-dbsize` (0: none).
    db_length: i64,
    /// The `DbSeqNum` option: `M` (0: none).
    db_num_seqs: i32,
    /// The volumes' search order (`RpsVolumeSet::order`), for `RpsVolumeSet::with_order`.
    order: Vec<usize>,
}

impl RpsUnits {
    /// The units of the database `db` with the options `-searchsp` (`searchsp`) and `-dbsize`
    /// (`dbsize`); `program` names a database error.
    pub(crate) fn of(db: &SeqDb, searchsp: i64, dbsize: i64, program: &str) -> Result<Self> {
        let volumes = RpsVolumeSet::of(db);
        let units = if volumes.order.len() <= 1 {
            vec![RpsUnit::of_name(db)]
        } else {
            volumes
                .order
                .iter()
                .map(|&index| RpsUnit::of_volume(&volumes.volumes[index]))
                .collect::<std::result::Result<Vec<_>, _>>()
                .map_err(|error| {
                    crate::algorithm::common::database::database_error(error, program)
                })?
        };
        let adjusted = (volumes.order.len() > 1 && searchsp == 0 && dbsize == 0).then(|| {
            let length = match db.total_length_stats() {
                0 => db.total_length(),
                stats => stats,
            };
            let count = match db.num_seqs_stats() {
                0 => db.ncbi_num_seqs(),
                stats => stats,
            };
            (length as i64, count)
        });
        let (db_length, db_num_seqs) = adjusted.unwrap_or((dbsize, 0));
        Ok(Self {
            units,
            adjusted,
            db_length,
            db_num_seqs,
            order: volumes.order,
        })
    }

    /// The volumes' search order (`RpsVolumeSet::order` of the database the units are of).
    pub(crate) fn order(&self) -> &[usize] {
        &self.order
    }

    /// The units in the order of their search.
    pub(crate) fn units(&self) -> &[RpsUnit] {
        &self.units
    }

    /// Whether the units are the volumes searched bare (a database of several volumes): no alias
    /// key reaches the search's Seq-aligns, which name each profile by its stored header
    /// (`report_ordinal_ids`).
    pub(crate) fn searches_volumes_bare(&self) -> bool {
        self.units.len() > 1
    }

    /// The first unit, whose statistics the merged result keeps and the reports print.
    pub(crate) fn first(&self) -> &RpsUnit {
        &self.units[0]
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-727 (cited at
    // `stats::setup::eff_lengths_database`)
    /// `BLAST_CalcEffLengths`' database length and number of sequences of `unit` (`T'`, `n'`):
    /// the `DbLength` option when above 0, else the setup's total; the `DbSeqNum` option when
    /// above 0, else the setup's count.
    pub(crate) fn eff_lengths(&self, unit: &RpsUnit) -> (i64, i32) {
        (
            if self.db_length > 0 {
                self.db_length
            } else {
                unit.setup_total
            },
            if self.db_num_seqs > 0 {
                self.db_num_seqs
            } else {
                unit.setup_count
            },
        )
    }

    /// The Gumbel block's database length of `unit` (blast_setup.c:912-923): the `DbLength` option
    /// when not 0 (also below 0), else the unit's total before the fallback to OID 0 (0 or
    /// negative for an alias file's `LENGTH` of 0 or of 2^63 or more: every Spouge e-value is
    /// then negative, written 0).
    pub(crate) fn gumbel_length(&self, unit: &RpsUnit) -> i64 {
        if self.db_length != 0 {
            self.db_length
        } else {
            unit.total
        }
    }

    /// The effective search space of a valid query context of `query_length` letters in
    /// `unit`'s search (`stats::setup::calc_eff_lengths` with `eff_lengths`, the Karlin block
    /// `ka`): `-searchsp` (`searchsp`) when not 0, and 0 where NCBI computes none (a length of 0).
    pub(crate) fn search_space(
        &self,
        unit: &RpsUnit,
        query_length: i64,
        ka: &crate::stats::KarlinParams,
        searchsp: i64,
        program: &str,
    ) -> Result<i64> {
        let (length, count) = self.eff_lengths(unit);
        if length == 0 && searchsp == 0 {
            return Ok(0);
        }
        Ok(crate::stats::setup::calc_eff_lengths(
            query_length,
            length,
            count,
            ka,
            searchsp,
            program,
        )?
        .eff_searchsp)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:155-181
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_memento_priv.hpp:107-147 (the
    // memento restores the local options after each unit's search; the remote options, which a
    // request holds, keep the value of the last unit)
    /// The `EffectiveSearchSpace` option that the last unit's search leaves in the options handle
    /// for a context of `query_length` letters (0 for an invalid context) without `-searchsp`
    /// (NCBI's request holds `-searchsp` then): `BLAST_CalcEffLengths` with the unit's total and
    /// count before the setup's fallback to OID 0 under the options (`DbLength`, `DbSeqNum`); 0
    /// for a length of 0, the context's query length for a negative one (no length adjustment,
    /// an effective database length of 1). The search strategy's requests after the first, and
    /// the BLAST archive's, hold it (`report::archive::OptionsHandle::pass`).
    pub(crate) fn request_space(
        &self,
        query_length: i64,
        ka: &crate::stats::KarlinParams,
        program: &str,
    ) -> Result<i64> {
        let unit = &self.units[self.units.len() - 1];
        match self.split_query_lengths(unit) {
            Some((length, count)) if query_length > 0 => Ok(crate::stats::setup::calc_eff_lengths(
                query_length,
                length,
                count,
                ka,
                0,
                program,
            )?
            .eff_searchsp),
            _ => Ok(0),
        }
    }

    /// The database length and count of `SplitQuery_SetEffectiveSearchSpace` for `unit`
    /// (cited above `request_space`): the unit's total before the setup's fallback to OID 0 and
    /// its count without the fallback, under the `DbLength` and `DbSeqNum` options; `None` where
    /// that length is 0 (`BLAST_CalcEffLengths` then sets no space, blast_setup.c:729-732).
    fn split_query_lengths(&self, unit: &RpsUnit) -> Option<(i64, i32)> {
        let length = if self.db_length > 0 {
            self.db_length
        } else {
            unit.total
        };
        let count = if self.db_num_seqs > 0 {
            self.db_num_seqs
        } else {
            unit.request_count
        };
        (length != 0).then_some((length, count))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-221 and
    // c++/src/algo/blast/core/blast_setup.c:775-778, 822-841
    /// The effective lengths of `unit`'s preliminary search where they differ from the setup's
    /// (`eff_lengths`), else `None`: NCBI's preliminary stage reads the search space that
    /// `SplitQuery_SetEffectiveSearchSpace` puts in the options for each context (under a memento
    /// that the traceback's set-up no longer sees), with `split_query_lengths`, and the setup's
    /// own where that is 0. Both give the same lengths where the unit's total is above 0 (and
    /// with `-searchsp` (`searchsp`), which both stages read), so they differ only over one
    /// volume of a total of 0 or less: a `LENGTH` of 2^63 or more (a negative total: the setup
    /// reads OID 0's length and 1, SplitQuery the negative total and the count), and with
    /// `-dbsize` a `LENGTH` of 0 too (SplitQuery `(dbsize, count)`, the setup `(dbsize, 1)`;
    /// review B of D1a-4, N2); the hit cutoff of the preliminary stage reads them without the
    /// Gumbel block (`RpsRun::prelim_spaces`).
    pub(crate) fn prelim_lengths(&self, unit: &RpsUnit, searchsp: i64) -> Option<(i64, i32)> {
        if searchsp != 0 {
            return None;
        }
        self.split_query_lengths(unit)
            .filter(|&lengths| lengths != self.eff_lengths(unit))
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:1002-1010 (cited at
    // `common::database::report_sum_refusal`): the serial writer's `-n` of INT64_MIN.
    /// The `DbLength` option that `x_AdjustDbSize` writes into the options handle at each query
    /// batch's set-up (`search_sets`): `D` for any value where it runs, `None` elsewhere. The
    /// requests after the first hold it (the first is written before the first batch's set-up);
    /// a `D` of exactly -2^63 (an alias file's `LENGTH` or `STATS_TOTLEN` of 2^63), which NCBI's
    /// serial writer prints as `-0` (a signed overflow), is written as its digits
    /// (PD-LOSAT-NCBI-DEFECTS approved exception 17, Int8 -2^63 written as digits: the Owner's
    /// answer to D1a-3's Q1, A, of 2026-11-17, DECISIONS OQ-1117 of line RG).
    pub(crate) fn strategy_db_length(&self) -> Option<i64> {
        self.adjusted.map(|(length, _)| length)
    }

    /// Whether the setup of a unit fails with NCBI's "The average subject length is too short"
    /// (`RpsUnit::setup_subject_length`, `gumbel`: whether the search has the Gumbel block).
    pub(crate) fn subject_length_too_short(&self, gumbel: bool) -> bool {
        self.units
            .iter()
            .any(|unit| unit.setup_subject_length(gumbel).is_none())
    }
}

impl RpsUnit {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980 (cited at
    // `stats::setup::setup_subject_length`)
    /// LOSAT's refusal of the unit where its setup counts 0 sequences of a positive total in a
    /// search without the Gumbel block (`gumbel` false: `-searchsp`, `OLD_FSC`): NCBI divides the
    /// total by that count (`min_subject_length`) and crashes (SIGFPE, DECISIONS 14 of line RG).
    /// Only the unit of one volume can (an alias file's `NSEQ 0`): a bare volume counts its
    /// profiles, and `DbSeqNum` does not reach the division (D1a4-7). The division is the one of
    /// `setup_subject_length`, which NCBI's setup of the unit makes only for a batch with a valid
    /// query (`CLocalBlast::Run`'s `CheckInternalData` before it): rpstblastn checks it at the
    /// first such batch that searches the unit (`search_window`). rpsblast refuses it once before
    /// any output where a query of the call is valid, by the setup's count, which is the first
    /// unit's, the one that can count 0 (`blastp::blast_engine::call_setup::setup_count_refused`;
    /// review C of D1a-4, MF1 and N2; review D, N2).
    pub(crate) fn check_setup_count(&self, gumbel: bool, program: &str) -> Result<()> {
        crate::algorithm::common::database::check_setup_count(
            self.setup_total,
            self.setup_count,
            gumbel,
            program,
        )
    }

    /// `BLAST_GapAlignSetUp`'s `min_subject_length` of the unit's setup
    /// (`stats::setup::setup_subject_length`: the unit's shortest length with the Gumbel block,
    /// `gumbel`, else its total over its count as NCBI's `Int4`), the one owner of both programs'
    /// check (review A of D1a-4, MF5); `None` where NCBI's setup fails with
    /// `stats::setup::SUBJECT_LENGTH_TOO_SHORT` (an alias file's `NSEQ` above the `LENGTH` over
    /// one volume), after its effective lengths and before the RPS engine runs
    /// (`check_engine_length`). A
    /// count of 0 of a positive total is LOSAT's refusal before (`check_setup_count`).
    pub(crate) fn setup_subject_length(&self, gumbel: bool) -> Option<i32> {
        crate::stats::setup::setup_subject_length(
            self.setup_total,
            self.setup_count,
            gumbel.then_some((self.min_length, false)),
        )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1129-1149 (cited at
    // `RpsUnit` and in `rps_search`'s module header) and c++/src/algo/blast/core/blast_extend.c:42-67
    /// LOSAT's refusals of the unit's concatenated length (`concat_length`), the one owner of both
    /// programs and of every unit, where NCBI's RPS engine (`s_RPSPreliminarySearchEngine`) meets
    /// them: at the unit's search of the first query batch that NCBI searches, after that
    /// batch's query set-up and the unit's setup (so an empty query file, or a first batch
    /// without letters, ends as NCBI's; review B of D1a-4, MF1). A `dbsize` above `INT4_MAX`:
    /// NCBI returns -3 with the profiles attached to its score block, which frees them twice
    /// (the run aborts). A length `Lc` whose diagonal table never stops growing:
    /// `diag_array_length` is an `Int4` that doubles from 1 while it is below `Lc` plus the word
    /// window (`window`): it ends at 2^30 where the sum is at most 2^30, else it doubles to 2^31
    /// (`INT_MIN`) and then to 0, and stays 0. A sum above `INT4_MAX` (an `Lc` within the window
    /// of it, such as `LENGTH` 2^31 - 4 with `NSEQ` 3) overflows NCBI's `Int4` addition, and the
    /// widened sum here is above 2^30. A bare volume of several meets them only with its own
    /// rows, one volume also through an alias file's `LENGTH` and `NSEQ`; an `Lc` below the rows,
    /// also 0 or negative, is searched as NCBI searches it (`RpsRun::concat_length`).
    pub(crate) fn check_engine_length(&self, window: i32, program: &str) -> Result<()> {
        if self.engine_size > i128::from(i32::MAX) {
            bail!(
                "an RPS-BLAST volume whose length and number of profiles (an alias file's LENGTH and NSEQ over one volume) sum to {}, above 2147483647 (NCBI BLAST+ 2.17.0 then frees the profiles twice and aborts), is not supported by LOSAT's {program}",
                self.engine_size
            );
        }
        let concatenated = self.concat_length();
        if i64::from(concatenated) + i64::from(window) > 1 << 30 {
            bail!(
                "an RPS-BLAST volume of {concatenated} concatenated positions (its length and number of profiles, or an alias file's LENGTH and NSEQ over one volume), whose diagonal table NCBI BLAST+ 2.17.0 never stops growing (above 2^30 entries with the window), is not supported by LOSAT's {program}"
            );
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:349-371,374-379
// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:73-83
/// What the set-up of each query batch of an RPS search (rpsblast, rpstblastn) writes into NCBI's
/// options handle, in its order, which the search strategy's requests after the first hold
/// (`report::archive::OptionsHandle::pass`; NCBI BLAST+ 2.17.0 `p9.rps.q40.f6`,
/// `p9.rps.3v.q40.n1.f6`, `p9.rtn.c0.bound1.f6`): a database of several volumes (names) without
/// `-searchsp` or `-dbsize` gets its length (`x_AdjustDbSize`, before the volumes; `DbSeqNum` is no
/// remote option): `db_length`, `RpsUnits::strategy_db_length`, then each volume's set-up the
/// matrix and the gap costs of its `.aux`
/// (`CreateRpsStructures`; the volumes' are the same, `RpsVolumeSetup`), also in a batch without a
/// valid query. `costs` are the volumes' matrix and gap costs as their files' set-up read them
/// (`RpsVolumeSetup::costs`), `None` when the files do not open: the batch's search fails there with
/// NCBI's error and no request follows. NCBI's threads over the volumes set them up on clones
/// (`RunThreadedSearch`), whose search spaces never reach the handle; LOSAT writes the requests of
/// NCBI's `-num_threads 1` at every thread count (DECISIONS F1b-2 of line RG).
pub(crate) fn search_sets(
    db_length: Option<i64>,
    costs: Option<(String, i32, i32)>,
) -> Vec<crate::report::archive::options::Param> {
    use crate::report::archive::fields;
    use crate::report::archive::options::{Param, ParamValue};
    let mut sets = Vec::new();
    if let Some(db_length) = db_length {
        // `SetDbLength` of `D` as an `Int8`, whatever its sign (0 too; probes p6_3v_l0#3,
        // p6_3v_l2p63p1#3 of line RG D1a-4).
        sets.push(Param::new(
            fields::DB_LENGTH,
            ParamValue::BigInteger(db_length),
        ));
    }
    if let Some((matrix, gap_open, gap_extend)) = costs {
        sets.push(Param::new(fields::MATRIX_NAME, ParamValue::Str(matrix)));
        sets.push(Param::new(
            fields::GAP_OPENING_COST,
            ParamValue::Integer(gap_open),
        ));
        sets.push(Param::new(
            fields::GAP_EXTENSION_COST,
            ParamValue::Integer(gap_extend),
        ));
    }
    sets
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1414-1418 (quoted
// at `SeqDb::included_ordinal`): the formatter's `CSeqDB`, opened over the `-db` name with its
// alias keys, finds no OID that the alias files' ranges leave out.
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:174-185 (cited in
// PD-LOSAT-NCBI-DEFECTS 1.37): `sacc` and `saccver` of such a subject are a null pointer exception.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:810-813 and
// c++/src/objtools/align_format/align_format_util.cpp:1014-1039 (`PruneSeqalign`)
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqinfosrc_aux.cpp:51-57 and
// c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:793-813 (the Seq-aligns of each volume's
// results name a profile by the Seq-ids of its header, read by OID whatever the ranges include)
/// The ordinal ids of the profiles of one query batch's final hit lists `lists`, renamed from
/// their OIDs (`shift_ordinal_ids`), for the reports, which look them up by ordinal id: NCBI's
/// lookup table holds every profile of a volume and its search ignores the alias files'
/// `FIRST_OID`, `LAST_OID` and `OIDLIST`, while its formatter does not find a profile that they leave out
/// (`Unknown`, "no longer exists in database", or a null pointer exception in `sacc` and
/// `saccver`). The formatter looks up only the profiles it reports: the first `cut` subjects of
/// each query's list (the hit list size, `PruneSeqalign` with `m_HitlistSize`, which counts
/// following profiles of one Seq-align id as one subject: `formatter_kept`; outfmt 0's
/// descriptions and alignments are within it), after outfmt 7 has counted the HSPs of the whole
/// list. A profile left out among them is LOSAT's refusal of the batch (the D1a-4 entry of
/// PD-LOSAT-NCBI-DEFECTS, the classes of 1.10 and 1.37: reported sequences; DECISIONS D1a4-4 of
/// line RG), one
/// text for `rpsblast` and `rpstblastn` (`program`). A profile left out past the cut (several
/// volumes, whose merged lists are longer than the hit list size) leaves its list: nothing reads
/// it again but its header, which NCBI's Seq-aligns read by OID (a header file that fails there
/// is the batch's error, as NCBI's), and its HSPs, which outfmt 7's count and the prefetch's
/// limit still count, and whose subject rows NCBI's prefetch reads (`Unlisted`; `None` where no
/// list lost a profile). Nothing changes where the ordinal ids are the OIDs
/// (`SeqDb::ordinal_ids_are_oids`: no range).
///
/// The alias files' `OID_MASK_TYPE` (line RG D1b-2) reaches the Seq-aligns of a database of one
/// volume, searched by the `-db` name: a profile it leaves without a line gets none. The volumes
/// of several are searched bare (`bare`, `RpsUnits::searches_volumes_bare`): their Seq-aligns
/// name every profile by its stored header, and the formatter, which reads the profiles over the
/// `-db` name, does not find a reported one whose lines the type removes (probe xp3v of D1b-2:
/// `Unknown`, "no longer exists in database", or a null pointer exception), LOSAT's refusal of
/// the batch in the same classes; one past the cut leaves its list as a profile out of the
/// ranges does.
pub(crate) fn report_ordinal_ids(
    db: &SeqDb,
    lists: &mut [Option<BlastpHitList>],
    cut: usize,
    program: &str,
    bare: bool,
) -> Result<Option<Unlisted>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1580-1586 (cited at
    // `SeqDb::drops_lineless_hits`)
    // In a database of one volume, searched by the `-db` name, a profile whose every definition
    // line the alias files' `OID_MASK_TYPE` erased gets no Seq-align (its header is read by OID,
    // whatever the ranges include): its lists leave the reports, outfmt 7's count and the
    // prefetch, before the formatter's cut (probes rps_* of line RG D1b-2: the hits of profile
    // `ref|XP_000101.1|` are gone, exit 0). NCBI reads the header of every listed profile there,
    // past the formatter's cut too: the first that cannot be read or that LOSAT refuses, in the
    // lists' order, is the batch's error. The volumes of several are searched bare (below; probe
    // xp3v of D1b-2: NCBI's Seq-aligns keep profile `ref|XP_000102.1|`).
    let erasing = db.drops_lineless_hits();
    if erasing && !bare {
        for list in lists.iter_mut().flatten() {
            list.hsplist_array.truncate(list.hsplist_count);
            let keep = list
                .hsplist_array
                .iter()
                .map(|hsp_list| db.keeps_line_of_oid(hsp_list.oid))
                .collect::<Result<Vec<bool>, _>>()
                .map_err(|error| {
                    crate::algorithm::common::database::database_error(error, program)
                })?;
            let mut keep = keep.into_iter();
            list.hsplist_array.retain(|_| keep.next().unwrap_or(true));
            list.hsplist_count = list.hsplist_array.len();
        }
    }
    if db.ordinal_ids_are_oids() && !(erasing && bare) {
        return Ok(None);
    }
    let queries = lists.len();
    let mut unlisted: Option<Unlisted> = None;
    // Whether a dropped profile's HSP has a subject row that starts with a gap.
    let starts_with_gap = |hsp_list: &crate::algorithm::blastp::hsp::BlastpHspList| {
        hsp_list.hsps.iter().any(|hsp| {
            super::redo::starts_with_row_gap(
                hsp,
                crate::report::serial::seqalign::AlignRow::Subject,
            )
        })
    };
    for (index, slot) in lists.iter_mut().enumerate() {
        let Some(list) = slot.as_mut() else {
            continue;
        };
        list.hsplist_array.truncate(list.hsplist_count);
        let cut = formatter_kept(db, list, cut, bare)
            .map_err(|error| crate::algorithm::common::database::database_error(error, program))?;
        let mut position = 0usize;
        let mut dropped_hsps = 0usize;
        let mut dropped_gap_start = false;
        let mut failed = None;
        list.hsplist_array.retain_mut(|hsp_list| {
            if failed.is_some() {
                return true;
            }
            let reported = position < cut;
            position += 1;
            // A profile of the bare volumes whose lines the alias files remove (their
            // `OID_MASK_TYPE` or their id lists, `SeqDb::drops_lineless_hits`): the formatter
            // does not find it among those it reports, and nothing reads it past the cut but its
            // stored header (the Seq-align's read; `erases_a_line_of_oid` reads the same bytes,
            // and meets NCBI's exception of an id list's lookup there).
            if erasing && bare {
                match db.erases_a_line_of_oid(hsp_list.oid) {
                    Err(error) => {
                        failed = Some(crate::algorithm::common::database::database_error(
                            error, program,
                        ));
                        return true;
                    }
                    Ok(true) if reported => {
                        failed = Some(anyhow::anyhow!(
                            "a hit of the profile of OID {} whose definition lines the RPS-BLAST database's alias files remove (their OID_MASK_TYPE, GILIST, TILIST, SEQIDLIST, TAXIDLIST or MEMB_BIT), in a database of several volumes (NCBI BLAST+ 2.17.0 searches the volumes without the alias files, and its report does not find the profile: Unknown, or a null pointer exception) is not supported by LOSAT's {program}",
                            hsp_list.oid
                        ));
                        return true;
                    }
                    Ok(true) => {
                        dropped_hsps += hsp_list.hsps.len();
                        dropped_gap_start |= starts_with_gap(hsp_list);
                        return false;
                    }
                    Ok(false) => {}
                }
            }
            // A reported profile out of the range is found where an included copy of its record
            // holds its Seq-align id (`found_as_copy`).
            let ordinal = match db.included_ordinal(hsp_list.oid) {
                None if reported => match found_as_copy(db, hsp_list.oid) {
                    Ok(Found::Copy(ordinal)) => Some(ordinal),
                    Ok(Found::Nothing) => {
                        failed = Some(anyhow::anyhow!(
                            "a hit of the profile of OID {} outside the RPS-BLAST database's alias files' FIRST_OID, LAST_OID and OIDLIST (NCBI BLAST+ 2.17.0's report does not find it: Unknown, or a null pointer exception) is not supported by LOSAT's {program}",
                            hsp_list.oid
                        ));
                        return true;
                    }
                    Ok(Found::Other(other)) => {
                        failed = Some(anyhow::anyhow!(
                            "a hit of the profile of OID {} outside the RPS-BLAST database's alias files' FIRST_OID, LAST_OID and OIDLIST, whose Seq-id NCBI BLAST+ 2.17.0's report finds in another record or other letters (the profile of OID {other}), is not supported by LOSAT's {program}",
                            hsp_list.oid
                        ));
                        return true;
                    }
                    Err(error) => {
                        failed = Some(crate::algorithm::common::database::database_error(
                            error, program,
                        ));
                        return true;
                    }
                },
                included => included,
            };
            match ordinal {
                Some(ordinal) => {
                    hsp_list.oid = ordinal;
                    for hsp in &mut hsp_list.hsps {
                        hsp.s_idx = ordinal;
                    }
                    true
                }
                None => {
                    let header = if bare {
                        db.read_stored_header_of_oid(hsp_list.oid)
                    } else {
                        db.read_header_of_oid(hsp_list.oid)
                    };
                    if let Err(error) = header {
                        failed = Some(crate::algorithm::common::database::database_error(
                            error, program,
                        ));
                        return true;
                    }
                    // Every HSP of the RPS traceback scores above 0: each is a Seq-align.
                    dropped_hsps += hsp_list.hsps.len();
                    dropped_gap_start |= starts_with_gap(hsp_list);
                    false
                }
            }
        });
        list.hsplist_count = list.hsplist_array.len();
        if let Some(error) = failed {
            return Err(error);
        }
        if dropped_hsps > 0 {
            let unlisted = unlisted.get_or_insert_with(|| Unlisted {
                hsps: vec![0; queries],
                subject_row_throws: false,
            });
            unlisted.hsps[index] = dropped_hsps;
            unlisted.subject_row_throws |= dropped_gap_start;
        }
    }
    Ok(unlisted)
}

/// The HSPs of the profiles that left the lists of a query batch past the formatter's cut
/// (`report_ordinal_ids`): NCBI's Seq-aligns still hold them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Unlisted {
    /// Their number in each query's list: outfmt 7's count and the prefetch's limit count them.
    pub(crate) hsps: Vec<usize>,
    /// Whether one of them starts with a gap in the subject's row
    /// (`serial::seqalign::gapped_row_range`): RPSTBLASTN's prefetch, which reads the subject row
    /// of every Std-seg Seq-align, throws there (`blastx::runtime::BatchResults::
    /// unlisted_subject_row_throws`). RPSBLAST's Seq-aligns are Dense-segs: it does not read it.
    pub(crate) subject_row_throws: bool,
}

/// What NCBI's reports find of a reported profile that the alias files' ranges leave out
/// (`found_as_copy`).
enum Found {
    /// An included copy of its record, under this ordinal id.
    Copy(u32),
    /// No included profile holds its Seq-align id: NCBI writes `Unknown` or stops with a null
    /// pointer exception.
    Nothing,
    /// An included profile of another record, or of other letters, holds its Seq-align id (this
    /// OID).
    Other(u32),
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1381-1459 (quoted
// at `seqdb::SeqidLookup::seqid_to_oids`): the formatter's lookup of a Seq-id gives the OIDs that
// the alias files' ranges include.
/// What NCBI's reports find of the profile of OID `oid`, which the alias files' ranges leave out:
/// the first included profile that holds its Seq-align id (`SeqidLookup::seqid_to_oids`), shown
/// in its place; LOSAT writes it where that profile is a copy of its record (the same definition
/// lines but for their taxids, the rule of `report::db_subject::joins_as_copy`, and the same
/// letters): a volume given twice under two paths, one of them past the range (probe dup6 of line
/// RG D1a-4 fix round 2: NCBI's rows of both copies, exit 0). The letters: NCBI's reports show the
/// found profile's length and letters (its scope's Bioseq for the Seq-align id) with this
/// profile's HSPs, whose identities LOSAT counts on this profile's consensus
/// (`traceback::formatter_identities`), so a profile of the same definition lines and other
/// letters is refused as another record is (`report::db_subject`'s `fetched` refuses such a
/// subject in blastp's reports; review C of D1a-4, N4).
fn found_as_copy(db: &SeqDb, oid: u32) -> Result<Found, crate::seqdb::SeqDbError> {
    use crate::seqdb::SequenceHeader;
    let SequenceHeader::Parsed(set) = db.read_header_of_oid(oid)? else {
        return Ok(Found::Nothing);
    };
    let Some(id) = set.alignment_id() else {
        return Ok(Found::Nothing);
    };
    let Some(&found) = db.seqid_lookup()?.seqid_to_oids(id, false)?.first() else {
        return Ok(Found::Nothing);
    };
    let record = |set: &crate::seqid::DeflineSet| -> Vec<crate::seqid::BlastDefline> {
        set.deflines
            .iter()
            .map(|defline| crate::seqid::BlastDefline {
                taxid: None,
                ..defline.clone()
            })
            .collect()
    };
    let letters = |oid: u32| -> Result<Vec<u8>, crate::seqdb::SeqDbError> {
        Ok(db
            .read_sequence_at(db.locate_oid(oid)?)?
            .protein(0)
            .to_vec())
    };
    let copy = match db.read_header_of_oid(found)? {
        SequenceHeader::Parsed(other) => {
            record(&set) == record(&other) && letters(oid)? == letters(found)?
        }
        SequenceHeader::Ordinal(_) => false,
    };
    Ok(match db.included_ordinal(found) {
        Some(ordinal) if copy => Found::Copy(ordinal),
        _ => Found::Other(found),
    })
}

/// The number of HSP lists at the start of `list` that NCBI's formatter keeps in its cut to `cut`
/// subjects (`PruneSeqalign`, cited at `report_ordinal_ids`): following lists whose profiles'
/// Seq-align ids match count as one subject, as NCBI's Seq-aligns name them (the headers read by
/// OID, `SeqDb::read_header_of_oid`). Two volumes that hold one profile (a volume given twice
/// under two paths) give such lists; the reports keep the same number
/// (`common::database::joined_cut`): both take the one count, `common::database::pruned_kept`.
/// A list of `cut` lists or fewer is kept whole, unread. The volumes searched bare (`bare`) name
/// their profiles by the stored headers (`SeqDb::read_stored_header_of_oid`).
fn formatter_kept(
    db: &SeqDb,
    list: &BlastpHitList,
    cut: usize,
    bare: bool,
) -> Result<usize, crate::seqdb::SeqDbError> {
    use crate::seqid::defline_set::{compare, Sic};
    let lists = &list.hsplist_array[..list.hsplist_count.min(list.hsplist_array.len())];
    if lists.len() <= cut {
        return Ok(lists.len());
    }
    crate::algorithm::common::database::pruned_kept(
        cut,
        lists.iter().map(|hsp_list| {
            let header = if bare {
                db.read_stored_header_of_oid(hsp_list.oid)?
            } else {
                db.read_header_of_oid(hsp_list.oid)?
            };
            Ok(match header {
                crate::seqdb::SequenceHeader::Parsed(set) => set.alignment_id().cloned(),
                crate::seqdb::SequenceHeader::Ordinal(_) => None,
            })
        }),
        |previous, id| {
            Ok(matches!((previous, id), (Some(a), Some(b)) if compare(a, b) == Sic::Yes))
        },
    )
}

/// The `.psq` letters of a volume's profiles (`TracebackProfile::letters`: the traceback's and
/// the redo's query, the reports' identities), read when a search first needs them: NCBI reads a
/// profile's sequence only when its list reaches the traceback, so a volume without `.psq`
/// searches a batch that has no list, or no valid query, as NCBI does (review B of PS12a, S3).
/// Each profile's letters are its own (`L_pin` of them, its `.pin` length), whatever its `.rps`
/// columns (series PR, item PS10e).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1176-1177
pub(crate) struct ProfileLetters<'a> {
    setup: &'a RpsVolumeSetup,
    volume: &'a Volume,
    /// The volume's first OID in the database.
    first_oid: u32,
    /// The volume's letters, and whether one of them is 28 or more (`profile::has_wide`).
    block: std::sync::OnceLock<(crate::seqdb::SequenceBlock, bool)>,
}

impl<'a> ProfileLetters<'a> {
    /// The letters of `volume`, whose first OID in the database is `first_oid`.
    pub(crate) fn new(setup: &'a RpsVolumeSetup, volume: &'a Volume, first_oid: u32) -> Self {
        Self {
            setup,
            volume,
            first_oid,
            block: std::sync::OnceLock::new(),
        }
    }

    /// The profile `oid` of the volume, whose offsets are `offsets`, as NCBI's traceback reads it
    /// (`TracebackProfile`).
    pub(crate) fn profile(
        &self,
        oid: u32,
        offsets: &[i32],
        program: &'static str,
    ) -> Result<super::profile::TracebackProfile<'_>> {
        let (block, wide) = self.block()?;
        Ok(super::profile::TracebackProfile::new(
            oid,
            self.first_oid,
            offsets,
            block.protein(oid as usize),
            *wide,
            program,
        ))
    }

    /// The letters of the profile `oid` of the volume.
    pub(crate) fn of(&self, oid: u32) -> Result<&[u8]> {
        Ok(self.block()?.0.protein(oid as usize))
    }

    /// The volume's letters, read once (threads that ask together may each read them; one copy
    /// stays), with whether one of them is 28 or more.
    fn block(&self) -> Result<&(crate::seqdb::SequenceBlock, bool)> {
        if let Some(block) = self.block.get() {
            return Ok(block);
        }
        let block = self
            .volume
            .read_block(0, self.volume.num_oids(), 0, 0)
            .map_err(|error| {
                crate::algorithm::common::database::database_error(error, self.setup.program)
            })?;
        let wide = super::profile::has_wide(block.protein_bytes());
        Ok(self.block.get_or_init(|| (block, wide)))
    }
}

/// The set-up of the volumes of one run: the program of LOSAT's refusals, the matrix and gap
/// costs of the first volume set up, and the volumes it regrouped.
#[derive(Debug, PartialEq)]
pub(crate) struct RpsVolumeSetup {
    program: &'static str,
    first_costs: Option<(String, i32, i32)>,
    /// The volumes (`canonical_volume_path`) whose cells `check` or `open` found not in groups of
    /// buckets in order and regrouped (`RpsVolumeFiles::is_regrouped`; series PR, item PS12m):
    /// `reopen` regroups them again without reading every hit. This set-up is the one owner of the
    /// fact: a run that reopens volumes keeps the set-up that checked them.
    regrouped: std::collections::BTreeSet<Vec<u8>>,
    /// The `-db` value of a database of one volume, which NCBI's error of a missing `.freq`
    /// names (`with_database`); `None`: each volume's path.
    database_name: Option<String>,
    /// The volumes (`canonical_volume_path`) whose files `check` or `open` has checked: `open`
    /// opens them again as `reopen` does (series PR, item PS13a: deltablast's domain search opens
    /// the volume for each query batch, with this set-up kept from search to search,
    /// `engine::RpsRunState`).
    checked: std::collections::BTreeSet<Vec<u8>>,
}

impl RpsVolumeSetup {
    pub(crate) fn new(program: &'static str) -> Self {
        Self {
            program,
            first_costs: None,
            regrouped: std::collections::BTreeSet::new(),
            database_name: None,
            checked: std::collections::BTreeSet::new(),
        }
    }

    /// The set-up of the volumes of the database `name` (the `-db` value as NCBI keeps it) of
    /// `num_volumes` volumes: NCBI searches a database of one volume by its name, so that the
    /// error of a missing `.freq` names it as given, and each volume of several by its path.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:341-346, 380-386
    pub(crate) fn with_database(mut self, name: &str, num_volumes: usize) -> Self {
        self.database_name = (num_volumes == 1).then(|| name.to_string());
        self
    }

    /// The files of `volume` (with its frequency ratios when `cbs`, composition-based
    /// statistics) after LOSAT's checks (`check`): `Ok(Err(error))` is NCBI's own error of the
    /// files, which the caller reports where NCBI's search raises it; `Err` is LOSAT's refusal.
    /// The profiles' letters come when the search needs them (`ProfileLetters`). A volume that this
    /// set-up has checked opens as `reopen` opens it: its files are those it checked (the run's),
    /// and NCBI's own errors of them come again where each search opens them.
    /// Uses `pool` for parallel `.loo` hit checking over the mapping, or reads the overflow array
    /// through a 1 MiB sliding buffer (`buffered_loo_check`) to keep unused `.loo` pages out of
    /// the resident set (`RpsVolumeFiles::open_for_in_pool`).
    pub(crate) fn open_in_pool(
        &mut self,
        volume: &Volume,
        cbs: bool,
        pool: Option<&crate::utils::threading::SearchPool<'_>>,
        buffered_loo_check: bool,
    ) -> Result<std::result::Result<RpsVolumeFiles, RpsVolumeError>> {
        let path = canonical_volume_path(volume.base());
        if self.checked.contains(&path) {
            return self.reopen(volume, cbs);
        }
        let files = match self.files_in_pool(volume, cbs, true, pool, buffered_loo_check)? {
            Ok(files) => files,
            Err(error) => return Ok(Err(error)),
        };
        self.check_files(volume, &files)?;
        self.checked.insert(path);
        Ok(Ok(files))
    }

    /// The matrix and gap costs of the volumes' `.aux` (`search_sets`): those of the first volume
    /// that `check` or `open` set up, else `volume`'s as its files read them
    /// (`RpsVolumeFiles::open_for`, without LOSAT's checks); `None` when its files do not open.
    pub(crate) fn costs(&self, volume: &Volume) -> Option<(String, i32, i32)> {
        if let Some(costs) = &self.first_costs {
            return Some(costs.clone());
        }
        match self.files(volume, false, true) {
            Ok(Ok(files)) => Some((
                files.aux.matrix.clone(),
                files.aux.gap_open,
                files.aux.gap_extend,
            )),
            _ => None,
        }
    }

    /// LOSAT's checks of `volume` alone, once per run before the first byte of output
    /// (rpstblastn, which maps the volumes again for each query batch: `reopen`).
    pub(crate) fn check(
        &mut self,
        volume: &Volume,
        cbs: bool,
    ) -> Result<std::result::Result<(), RpsVolumeError>> {
        let files = match self.files(volume, cbs, true)? {
            Ok(files) => files,
            Err(error) => return Ok(Err(error)),
        };
        self.check_files(volume, &files)?;
        self.checked.insert(canonical_volume_path(volume.base()));
        Ok(Ok(()))
    }

    /// The files of `volume`, which `check` has checked in this run, opened again for a search
    /// (the profiles' letters: `ProfileLetters`) without LOSAT's checks of `check` (series PR, item
    /// PS12s, D6: `RpsVolumeFiles::open_checked`); NCBI's own errors of the files come where its
    /// search of a batch meets them.
    pub(crate) fn reopen(
        &self,
        volume: &Volume,
        cbs: bool,
    ) -> Result<std::result::Result<RpsVolumeFiles, RpsVolumeError>> {
        self.files(volume, cbs, false)
    }

    /// The volume's files (and frequency ratios when `cbs`), the table's hits checked when
    /// `check_hits` (`RpsVolumeFiles::open_for`, `open_checked`).
    fn files(
        &self,
        volume: &Volume,
        cbs: bool,
        check_hits: bool,
    ) -> Result<std::result::Result<RpsVolumeFiles, RpsVolumeError>> {
        self.files_in_pool(volume, cbs, check_hits, None, false)
    }

    fn files_in_pool(
        &self,
        volume: &Volume,
        cbs: bool,
        check_hits: bool,
        pool: Option<&crate::utils::threading::SearchPool<'_>>,
        buffered_loo_check: bool,
    ) -> Result<std::result::Result<RpsVolumeFiles, RpsVolumeError>> {
        let num_seqs = volume.num_oids();
        let lengths = (0..num_seqs)
            .map(|oid| volume.protein_length(oid))
            .collect::<std::result::Result<Vec<u32>, _>>()
            .map_err(|error| {
                crate::algorithm::common::database::database_error(error, self.program)
            })?;
        let opened = if check_hits {
            RpsVolumeFiles::open_for_in_pool(
                volume.base(),
                &lengths,
                self.program,
                pool,
                buffered_loo_check,
            )
        } else {
            let regroup = self
                .regrouped
                .contains(&canonical_volume_path(volume.base()));
            RpsVolumeFiles::open_checked(volume.base(), &lengths, self.program, regroup)
        };
        let mut files = match opened {
            Ok(files) => files,
            Err(RpsVolumeFault::Files(RpsOpenError::Refused(text))) => bail!(text),
            Err(RpsVolumeFault::Files(RpsOpenError::Ncbi(text))) => {
                return Ok(Err(RpsVolumeError::Files(text)))
            }
            Err(RpsVolumeFault::NoMatrix) => return Ok(Err(RpsVolumeError::NoMatrix)),
        };
        if cbs {
            match files.open_freq(volume.base(), self.program, self.database_name.as_deref()) {
                Ok(()) => {}
                Err(RpsOpenError::Refused(text)) => bail!(text),
                Err(RpsOpenError::Ncbi(text)) => return Ok(Err(RpsVolumeError::Files(text))),
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:72-84
        // The search's Karlin, Gumbel and cutoff values are those of BLOSUM62 11/1 (blastp's
        // task); a database of other costs is refused, and so are volumes whose costs differ
        // (each volume sets its own in NCBI; the statistics are one set).
        Ok(Ok(files))
    }

    /// LOSAT's refusals of what NCBI's search does not reproducibly run on the volume's files.
    fn check_files(&mut self, volume: &Volume, files: &RpsVolumeFiles) -> Result<()> {
        if files.is_regrouped() {
            self.regrouped.insert(canonical_volume_path(volume.base()));
        }
        let aux = &files.aux;
        if aux.matrix != "BLOSUM62" || aux.gap_open != 11 || aux.gap_extend != 1 {
            bail!(
                "a profile database of the matrix {} with gap costs {}/{} is not supported by LOSAT's {}",
                aux.matrix,
                aux.gap_open,
                aux.gap_extend,
                self.program
            );
        }
        let costs = (aux.matrix.clone(), aux.gap_open, aux.gap_extend);
        if self.first_costs.get_or_insert_with(|| costs.clone()) != &costs {
            bail!(
                "a profile database whose volumes differ in the matrix or the gap costs is not supported by LOSAT's {}",
                self.program
            );
        }
        // The scores of an alignment of a profile add up to at most the sum of its rows' best
        // scores; above `INT4_MAX` NCBI's 32-bit scores overflow (undefined behavior). The
        // file is read through a small buffer, not the mapping.
        let bounds = files.profile_score_bounds(volume.base()).ok_or_else(|| {
            anyhow::anyhow!(
                "a profile database whose .rps file cannot be read is not supported by LOSAT's {}",
                self.program
            )
        })?;
        if let Some(profile) = bounds.iter().position(|&bound| bound > i64::from(i32::MAX)) {
            bail!(
                "a profile database whose profile {profile} can score {} (above 2147483647, where NCBI BLAST+'s 32-bit scores overflow) is not supported by LOSAT's {}",
                bounds[profile],
                self.program
            );
        }
        // NCBI's scaled values: `scale` cast to `Int4` multiplies the costs and cutoffs (a whole
        // scale of 1 or more here; `RpsRun::search_volume` refuses one whose products leave `Int4`
        // and one that changes Lambda).
        let scale = aux.scale;
        if !(scale >= 1.0 && scale.fract() == 0.0) {
            bail!(
                "a profile database of the scale {scale} (not a whole number of 1 or more) is not supported by LOSAT's {}",
                self.program
            );
        }
        // (The concatenated length is refused, or not, by its unit where NCBI's engine first
        // meets it: `RpsUnit::check_engine_length`. The profiles' letters are read as stored,
        // `Blast_TracebackGetEncoding` gives `eBlastEncodingProtein`: `ProfileLetters`.)
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:177-190
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1336-1346
// The lookup table is made after the score block, so the queries' Karlin-Altschul values
// keep their letters; the search (scan, extensions, traceback, rescale) reads them masked.
pub(crate) fn mask_unsupported_letters(aa_seq: &[u8], alphabet: usize) -> Option<Vec<u8>> {
    const PROT_MASK: u8 = 21;
    let inner = &aa_seq[1..aa_seq.len() - 1];
    (alphabet < crate::core::blast_util::BLASTAA_SIZE
        && inner.iter().any(|&letter| usize::from(letter) >= alphabet))
    .then(|| {
        let mut copy = aa_seq.to_vec();
        let end = copy.len() - 1;
        for letter in &mut copy[1..end] {
            if usize::from(*letter) >= alphabet {
                *letter = PROT_MASK;
            }
        }
        copy
    })
}

/// A volume's local ordinal ids as the database's: `start` is the volume's first ordinal id.
pub(crate) fn shift_ordinal_ids(hit_list: &mut BlastpHitList, start: u32) {
    if start == 0 {
        return;
    }
    let count = hit_list.hsplist_count;
    for list in hit_list.hsplist_array.iter_mut().take(count) {
        list.oid += start;
        for hsp in &mut list.hsps {
            hsp.s_idx += start;
        }
    }
}

/// The volumes after the first merge into the first's hit lists, query by query and volume by
/// volume, in the order of `s_MergeAlignSet`: the lists are the groups of Seq-aligns of one
/// subject (in the order of the Seq-aligns, `hit_list_to_seqalign_order`); the first Seq-align of
/// the current group of each side decides, the smaller E-value first, at equal E-values the
/// larger bit score, at equal both the final side's; a group of the input goes in whole before
/// the final group it beats, and when the final side runs out the rest of the input is
/// appended. Nothing is cut: the formatter keeps the first `max_target_seqs` subjects of the
/// merged set, and outfmt 7 counts its HSPs before the cut.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:69-137
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:155-166
#[cold]
#[inline(never)]
pub(crate) fn merge_volume(
    hit_lists: &mut [Option<BlastpHitList>],
    volume: Vec<Option<BlastpHitList>>,
) {
    for (slot, input) in hit_lists.iter_mut().zip(volume) {
        let Some(mut input) = input else {
            continue;
        };
        crate::algorithm::blastp::hsp::hit_list_to_seqalign_order(&mut input);
        if input.hsplist_count == 0 {
            continue;
        }
        match slot.as_mut() {
            None => {
                input.hsplist_array.shrink_to_fit();
                *slot = Some(input);
            }
            Some(final_list) => {
                crate::algorithm::blastp::hsp::hit_list_to_seqalign_order(final_list);
                if final_list.hsplist_count == 0 {
                    input.hsplist_array.shrink_to_fit();
                    *slot = Some(input);
                    continue;
                }
                merge_align_set(final_list, input);
            }
        }
    }
}

/// The E-value of a Seq-align's score list: 0 below `SMALLEST_EVALUE` (what `s_MergeAlignSet`
/// reads with `GetNamedScore`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1131 (and `:62`)
fn seqalign_evalue(e_value: f64) -> f64 {
    if e_value < 1.0e-180 {
        0.0
    } else {
        e_value
    }
}

/// `s_MergeAlignSet` of one query's lists (see [`merge_volume`]).
fn merge_align_set(final_list: &mut BlastpHitList, mut input: BlastpHitList) {
    final_list.hsplist_array.truncate(final_list.hsplist_count);
    input.hsplist_array.truncate(input.hsplist_count);
    let mut incoming = input.hsplist_array.into_iter().peekable();
    let mut at = 0usize;
    while let Some(group) = incoming.peek() {
        let (final_first, input_first) = (&final_list.hsplist_array[at].hsps[0], &group.hsps[0]);
        // The Seq-align carries 0 for an E-value below `SMALLEST_EVALUE`, so those compare equal.
        let (mut final_value, mut input_value) = (
            seqalign_evalue(final_first.e_value),
            seqalign_evalue(input_first.e_value),
        );
        if input_value == final_value {
            // The trick of NCBI's loop: the bit scores take the e-values' places.
            input_value = final_first.bit_score;
            final_value = input_first.bit_score;
        }
        if input_value < final_value {
            let group = incoming.next().expect("peeked");
            final_list.hsplist_array.insert(at, group);
            at += 1;
        } else {
            at += 1;
            if at == final_list.hsplist_array.len() {
                final_list.hsplist_array.extend(incoming);
                break;
            }
        }
    }
    final_list.hsplist_count = final_list.hsplist_array.len();
    final_list.hsplist_array.shrink_to_fit();
}

#[cfg(test)]
#[path = "volumes_tests.rs"]
mod tests;
