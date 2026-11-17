//! The volumes of a profile database as `rpsblast` and `rpstblastn` search them (series PR,
//! items PS11b and PS12a): NCBI's order of the volumes and their first ordinal ids, the search
//! units and their statistics under the alias files' keys (`RpsUnits`, line RG D1a-4), the set-up
//! of one volume (its files, LOSAT's refusals of what NCBI's search does not reproducibly run, the
//! consensus letters of its profiles), the masking of the letters that a profile has no column
//! for, the merge of the volumes' results and their ordinal ids for the reports. One owner for
//! both programs (review A of PS12a, N3).

use anyhow::{bail, Result};

use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::seqdb::rps::{RpsOpenError, RpsVolumeFault, RpsVolumeFiles};

/// NCBI's own errors of a volume's files, which its search of a query batch raises (`BLAST
/// engine error`, exit status 3).
pub(crate) enum RpsVolumeError {
    /// NCBI's text of a missing, short or unreadable file.
    Files(String),
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:905-937 (quoted at
// `SeqDb::statistics_total_length`) and c++/src/algo/blast/core/blast_engine.c:1129-1143
// ```c
//     num_db_seqs = BlastSeqSrcGetNumSeqs(seq_src);
//     dbsize = BlastSeqSrcGetTotLen(seq_src) + num_db_seqs;
//     if (dbsize > INT4_MAX)
//       return -3;
//     ...
//     concat_db.length = (Int4) dbsize;
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:131-139 (a bare
// volume walks no alias file: its index files store no shortest length, and it has no STATS key)
// ```c++
//             if ( m_MinLength <= 0)  m_MinLength = BLAST_SEQSRC_MINLENGTH;
//         }
//         m_NumSeqsStats     = x_GetNumSeqsStats();
//         m_TotalLengthStats = x_GetTotalLengthStats();
// ```
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
    /// The RPS engine's `dbsize`, `BlastSeqSrcGetTotLen + BlastSeqSrcGetNumSeqs` (NCBI's `Int8`,
    /// widened: it cannot overflow here), whose `Int4` is the concatenated length `Lc`
    /// (`concat_length`).
    pub(crate) engine_size: i128,
    /// Whether the unit's `BlastSeqSrcGetTotLen` (`GetTotalLength`, no STATS key) is above 0, the
    /// gate of the sort after the traceback's filters (`RpsRun::sorts_results`): not for an alias
    /// file's `LENGTH` of 0 or of 2^63 or more over one volume.
    pub(crate) sorts_results: bool,
}

impl RpsUnit {
    /// The concatenated length `Lc` of the unit's search, `concat_db.length = (Int4) dbsize`
    /// (blast_engine.c:1143, quoted above): the low 32 bits of `engine_size`. A bare volume's is
    /// the letters of its index and one delimiter row per profile, the rows of its `.rps` file; an
    /// alias file's `LENGTH` and `NSEQ` over one volume move it below the rows or above them, also
    /// to a negative value. The RPS engine reads it in three places (`RpsRun::concat_length`): the
    /// size of the diagonal table, the bound of the ungapped extensions to the right and the root
    /// of the gapped stage's interval tree.
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
            (i64::from(volume.protein_length(0)?), 1)
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
// ```c++
// CRef<CSearchResultSet> s_RunLocalRpsSearch(const string & db,
// 										   CBlastQueryVector  & query_vector,
// 										   CRef<CBlastOptionsHandle> opt_handle)
// {
// 	CSearchDatabase			search_db(db, CSearchDatabase::eBlastDbIsProtein);
// 	CRef<CLocalDbAdapter> 	db_adapter(new CLocalDbAdapter(search_db));
// ...
// 		if(1 == m_num_of_dbs)
// 		{
// 			return s_RunLocalRpsSearch(m_db_name, *m_query_vector, m_opt_handle);
// 		}
// 		else
// 		{
// 		   	s_ModifyVolumePaths(m_rps_databases);
// 			...
// 				results.push_back(s_RunLocalRpsSearch(m_rps_databases[i],
// 													  *m_query_vector,
// 													   m_opt_handle));
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:349-371 (quoted above
// `search_sets`): `x_AdjustDbSize`, before the volumes of several.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-727 (quoted at
// `stats::setup::eff_lengths_database`) and 912-923 (the Gumbel block's `db_length`, quoted at
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

    /// The first unit, whose statistics the merged result keeps and the reports print.
    pub(crate) fn first(&self) -> &RpsUnit {
        &self.units[0]
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-727 (quoted at
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
    // ```c++
    //     if (options->GetEffectiveSearchSpace() != 0) {
    //         return;
    //     }
    //
    //     const BlastSeqSrc* seqsrc = full_data->m_SeqSrc->GetPointer();
    //     Int8 total_length = BlastSeqSrcGetTotLenStats(seqsrc);
    //     if (total_length <= 0)
    //         total_length = BlastSeqSrcGetTotLen(seqsrc);
    //     Int4 num_seqs = BlastSeqSrcGetNumSeqsStats(seqsrc);
    //     if (num_seqs <= 0)
    //         num_seqs = BlastSeqSrcGetNumSeqs(seqsrc);
    //
    //     CEffectiveSearchSpaceCalculator calc(full_query_fact, *options,
    //                                          num_seqs, total_length,
    //                                          full_data->m_ScoreBlk->GetPointer());
    //     ...
    //     options->SetEffectiveSearchSpace(eff_searchsp);
    // ```
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
    /// (quoted above `request_space`): the unit's total before the setup's fallback to OID 0 and
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
    // ```c++
    //     CEffectiveSearchSpacesMemento eff_memento(m_Options);
    //     SplitQuery_SetEffectiveSearchSpace(m_Options, m_QueryFactory,
    //                                        m_InternalData);
    // ...
    //       Int8 effective_search_space =
    //           s_GetEffectiveSearchSpaceForContext(eff_len_options, index,
    //                                               blast_message);
    // ...
    //          if (effective_search_space == 0) {
    // ```
    /// The effective lengths of `unit`'s preliminary search where they differ from the setup's
    /// (`eff_lengths`), else `None`: NCBI's preliminary stage reads the search space that
    /// `SplitQuery_SetEffectiveSearchSpace` puts in the options for each context (under a memento
    /// that the traceback's set-up no longer sees), with `split_query_lengths`, and the setup's
    /// own where that is 0. Both give the same lengths where the unit's total is above 0 (and
    /// with `-searchsp` (`searchsp`), which both stages read), so they differ only for one
    /// volume whose `LENGTH` is 2^63 or more (a negative total: the setup reads OID 0's length
    /// and 1, SplitQuery the negative total and the count); the hit cutoff of the preliminary
    /// stage reads them without the Gumbel block (`RpsRun::prelim_spaces`).
    pub(crate) fn prelim_lengths(&self, unit: &RpsUnit, searchsp: i64) -> Option<(i64, i32)> {
        if searchsp != 0 {
            return None;
        }
        self.split_query_lengths(unit)
            .filter(|&lengths| lengths != self.eff_lengths(unit))
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:1002-1010 (quoted at
    // `common::database::report_sum_refusal`): the serial writer's `-n` of INT64_MIN.
    /// The `DbLength` option that `x_AdjustDbSize` writes into the options handle at each query
    /// batch's set-up (`search_sets`): `D` for any value where it runs, `None` elsewhere. A `D` of
    /// exactly -2^63 (an alias file's `LENGTH` or `STATS_TOTLEN` of 2^63), which NCBI's serial
    /// writer prints as `-0` (a signed overflow), is LOSAT's refusal, with the reason of the XML
    /// reports' (`common::database::STRATEGY_MINUS_ZERO_LENGTH`; DECISIONS D1a3-2 and D1a4-6 of
    /// line RG: until the Owner answers D1a-3's Q1).
    pub(crate) fn strategy_db_length(&self, program: &str) -> Result<Option<i64>> {
        match self.adjusted {
            Some((i64::MIN, _)) => Err(crate::algorithm::common::database::refusal(
                crate::algorithm::common::database::STRATEGY_MINUS_ZERO_LENGTH,
                program,
            )),
            adjusted => Ok(adjusted.map(|(length, _)| length)),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980 (quoted at
    // `stats::setup::setup_subject_length`)
    /// LOSAT's refusal of a unit whose setup counts 0 sequences of a positive total in a search
    /// without the Gumbel block (`gumbel` false: `-searchsp`, `OLD_FSC`): NCBI divides the total by
    /// that count (`min_subject_length`) and crashes (SIGFPE, DECISIONS 14 of line RG). Only the
    /// unit of one volume can (an alias file's `NSEQ 0`): a bare volume counts its profiles, and
    /// `DbSeqNum` does not reach the division (D1a4-7).
    pub(crate) fn check_setup_counts(&self, gumbel: bool, program: &str) -> Result<()> {
        if !gumbel
            && self
                .units
                .iter()
                .any(|unit| unit.setup_total > 0 && unit.setup_count == 0)
        {
            return Err(crate::algorithm::common::database::refusal(
                crate::algorithm::common::database::STATISTICS_COUNT_ZERO,
                program,
            ));
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1129-1143 (quoted at
    // `RpsUnit`) and c++/src/algo/blast/core/blast_extend.c:42-67
    // ```c
    // static BLAST_DiagTable*
    // s_BlastDiagTableNew (Int4 qlen, Boolean multiple_hits, Int4 window_size)
    // {
    //         ...
    //         diag_array_length = 1;
    //         /* What power of 2 is just longer than the query? */
    //         while (diag_array_length < (qlen+window_size)) {
    //                 diag_array_length = diag_array_length << 1;
    //         }
    // ```
    /// LOSAT's refusals of the concatenated length of a unit (`RpsUnit::concat_length`), the one
    /// owner of both programs and of every unit, checked once per run before the first byte of
    /// output: a `dbsize` above `INT4_MAX` (NCBI returns -3 with the profiles attached to its score
    /// block, which frees them twice: the run aborts), and a length `Lc` whose diagonal table never
    /// stops growing. `diag_array_length` is an `Int4` that doubles from 1 while it is below `Lc`
    /// plus the word window (`window`): it ends at 2^30 where the sum is at most 2^30, else it
    /// doubles to 2^31 (`INT_MIN`) and then to 0, and stays 0. A sum above `INT4_MAX` (an `Lc`
    /// within the window of it, such as `LENGTH` 2^31 - 4 with `NSEQ` 3) overflows NCBI's `Int4`
    /// addition, and the widened sum here is above 2^30. A bare volume of several meets them only
    /// with its own rows, one volume also through an alias file's `LENGTH` and `NSEQ`; an `Lc`
    /// below the rows, also 0 or negative, is searched as NCBI searches it
    /// (`RpsRun::concat_length`).
    pub(crate) fn check_engine_lengths(&self, window: i32, program: &str) -> Result<()> {
        for unit in &self.units {
            if unit.engine_size > i128::from(i32::MAX) {
                bail!(
                    "an RPS-BLAST volume whose length and number of profiles (an alias file's LENGTH and NSEQ over one volume) sum to {}, above 2147483647 (NCBI BLAST+ 2.17.0 then frees the profiles twice and aborts), is not supported by LOSAT's {program}",
                    unit.engine_size
                );
            }
            let concatenated = unit.concat_length();
            if i64::from(concatenated) + i64::from(window) > 1 << 30 {
                bail!(
                    "an RPS-BLAST volume of {concatenated} concatenated positions (its length and number of profiles, or an alias file's LENGTH and NSEQ over one volume), whose diagonal table NCBI BLAST+ 2.17.0 never stops growing (above 2^30 entries with the window), is not supported by LOSAT's {program}"
                );
            }
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:349-371,374-379
// ```c++
// void CLocalRPSBlast::x_AdjustDbSize(void)
// {
// 	if(m_opt_handle->GetOptions().GetEffectiveSearchSpace()!= 0)
// 		return;
//
// 	if(m_opt_handle->GetOptions().GetDbLength()!= 0)
// 		return;
//
// 	CSeqDB db(m_db_name, CSeqDB::eProtein);
//
// 	Uint8 db_size = db.GetTotalLengthStats();
// 	...
// 	if(0 == db_size)
// 	    db_size = db.GetTotalLength();
// 	...
// 	m_opt_handle->SetOptions().SetDbLength(db_size);
// 	m_opt_handle->SetOptions().SetDbSeqNum(num_seq);
// ...
// CRef<CSearchResultSet> CLocalRPSBlast::Run(void)
// {
// 	if(1 != m_num_of_dbs)
// 	{
// 		x_AdjustDbSize();
// 	}
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:73-83
// ```c++
// CSetupFactory::CreateRpsStructures(const string& rps_dbname,
//                                    CRef<CBlastOptions> options)
// {
//     ...
//     CRef<CBlastRPSInfo> retval(new CBlastRPSInfo(rps_dbname, mode));
//     options->SetMatrixName(retval->GetMatrixName());
//     options->SetGapOpeningCost(retval->GetGapOpeningCost());
//     options->SetGapExtensionCost(retval->GetGapExtensionCost());
//     return retval;
// }
// ```
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
    use crate::report::archive::options::{Param, ParamValue};
    let mut sets = Vec::new();
    if let Some(db_length) = db_length {
        // `SetDbLength(Int8)` of `D`, whatever its sign (0 too; probes p6_3v_l0#3,
        // p6_3v_l2p63p1#3 of line RG D1a-4).
        sets.push(Param {
            name: "DbLength",
            value: ParamValue::BigInteger(db_length),
        });
    }
    if let Some((matrix, gap_open, gap_extend)) = costs {
        sets.push(Param {
            name: "MatrixName",
            value: ParamValue::Str(matrix),
        });
        sets.push(Param {
            name: "GapOpeningCost",
            value: ParamValue::Integer(gap_open),
        });
        sets.push(Param {
            name: "GapExtensionCost",
            value: ParamValue::Integer(gap_extend),
        });
    }
    sets
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1414-1418 (quoted
// at `SeqDb::included_ordinal`): the formatter's `CSeqDB`, opened over the `-db` name with its
// alias keys, finds no OID that the alias files' ranges leave out.
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:174-185 (quoted in
// PD-LOSAT-NCBI-DEFECTS 1.37): `sacc` and `saccver` of such a subject are a null pointer exception.
/// The ordinal ids of the profiles of one query batch's final hit lists `lists`, renamed from
/// their OIDs (`shift_ordinal_ids`), for the reports, which look them up by ordinal id: NCBI's
/// lookup table holds every profile of a volume and its search ignores the alias files'
/// `FIRST_OID` and `LAST_OID`, while its formatter does not find a profile that they leave out
/// (`Unknown`, "no longer exists in database", or a null pointer exception in `sacc` and
/// `saccver`). Such a profile in the lists is LOSAT's refusal of the batch (PD-LOSAT-NCBI-DEFECTS
/// 1.40, the classes of 1.10 and 1.37; DECISIONS D1a4-4 of line RG), one text for `rpsblast` and
/// `rpstblastn` (`program`). Nothing changes where the ordinal ids are the OIDs
/// (`SeqDb::ordinal_ids_are_oids`: no range).
pub(crate) fn report_ordinal_ids(
    db: &SeqDb,
    lists: &mut [Option<BlastpHitList>],
    program: &str,
) -> Result<()> {
    if db.ordinal_ids_are_oids() {
        return Ok(());
    }
    for list in lists.iter_mut().flatten() {
        let count = list.hsplist_count.min(list.hsplist_array.len());
        for hsp_list in &mut list.hsplist_array[..count] {
            let Some(ordinal) = db.included_ordinal(hsp_list.oid) else {
                bail!(
                    "a hit of the profile of OID {} outside the RPS-BLAST database's alias file FIRST_OID and LAST_OID (NCBI BLAST+ 2.17.0's report does not find it: Unknown, or a null pointer exception) is not supported by LOSAT's {program}",
                    hsp_list.oid
                );
            };
            hsp_list.oid = ordinal;
            for hsp in &mut hsp_list.hsps {
                hsp.s_idx = ordinal;
            }
        }
    }
    Ok(())
}

/// The consensus letters of a volume's profiles at their offsets (the traceback's letters and the
/// reports' identities), read when a search first needs them: NCBI reads a profile's sequence
/// only when its list reaches the traceback, so a volume without `.psq` searches a batch that
/// has no list, or no valid query, as NCBI does (review B of PS12a, S3).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1176-1177
/// ```c
///       if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0)
///           continue;
/// ```
pub(crate) struct Consensus<'a> {
    setup: &'a RpsVolumeSetup,
    volume: &'a Volume,
    files: &'a RpsVolumeFiles,
    letters: std::sync::OnceLock<Vec<u8>>,
}

impl<'a> Consensus<'a> {
    pub(crate) fn new(
        setup: &'a RpsVolumeSetup,
        volume: &'a Volume,
        files: &'a RpsVolumeFiles,
    ) -> Self {
        Self {
            setup,
            volume,
            files,
            letters: std::sync::OnceLock::new(),
        }
    }

    /// The letters, read once (threads that ask together may each read them; one copy stays).
    pub(crate) fn letters(&self) -> Result<&[u8]> {
        if let Some(letters) = self.letters.get() {
            return Ok(letters);
        }
        let letters = self.setup.consensus(self.volume, self.files)?;
        Ok(self.letters.get_or_init(|| letters))
    }
}

/// The set-up of the volumes of one run: the program of LOSAT's refusals and the matrix and gap
/// costs of the first volume set up.
pub(crate) struct RpsVolumeSetup {
    program: &'static str,
    first_costs: Option<(String, i32, i32)>,
}

impl RpsVolumeSetup {
    pub(crate) fn new(program: &'static str) -> Self {
        Self {
            program,
            first_costs: None,
        }
    }

    /// The files of `volume` (with its frequency ratios when `cbs`, composition-based
    /// statistics) after LOSAT's checks (`check`): `Ok(Err(error))` is NCBI's own error of the
    /// files, which the caller reports where NCBI's search raises it; `Err` is LOSAT's refusal.
    /// The consensus letters come when the search needs them (`Consensus`).
    pub(crate) fn open(
        &mut self,
        volume: &Volume,
        cbs: bool,
    ) -> Result<std::result::Result<RpsVolumeFiles, RpsVolumeError>> {
        let files = match self.files(volume, cbs, true)? {
            Ok(files) => files,
            Err(error) => return Ok(Err(error)),
        };
        self.check_files(volume, &files)?;
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
        Ok(Ok(()))
    }

    /// The files of `volume`, which `check` has checked in this run, opened again for a search
    /// (the consensus letters: `Consensus`) without LOSAT's checks of `check` (series PR, item
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
        let num_seqs = volume.num_oids();
        let lengths = (0..num_seqs)
            .map(|oid| volume.protein_length(oid))
            .collect::<std::result::Result<Vec<u32>, _>>()
            .map_err(|error| {
                crate::algorithm::common::database::database_error(error, self.program)
            })?;
        let opened = if check_hits {
            RpsVolumeFiles::open_for(volume.base(), &lengths, self.program)
        } else {
            RpsVolumeFiles::open_checked(volume.base(), &lengths, self.program)
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
            match files.open_freq(volume.base()) {
                Ok(()) => {}
                Err(RpsOpenError::Refused(text)) => bail!(text),
                Err(RpsOpenError::Ncbi(text)) => return Ok(Err(RpsVolumeError::Files(text))),
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:72-83
        // ```c++
        //     if (rps_info) {
        //         opts_memento->m_ScoringOpts->matrix = strdup(rps_info->aux_info.orig_score_matrix);
        //         opts_memento->m_ScoringOpts->gap_open = rps_info->aux_info.gap_open_penalty;
        //         opts_memento->m_ScoringOpts->gap_extend = rps_info->aux_info.gap_extend_penalty;
        // ```
        // The search's Karlin, Gumbel and cutoff values are those of BLOSUM62 11/1 (blastp's
        // task); a database of other costs is refused, and so are volumes whose costs differ
        // (each volume sets its own in NCBI; the statistics are one set).
        Ok(Ok(files))
    }

    /// LOSAT's refusals of what NCBI's search does not reproducibly run on the volume's files.
    fn check_files(&mut self, volume: &Volume, files: &RpsVolumeFiles) -> Result<()> {
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
        // NCBI's scaled values: `(Int4)scale` multiplies the costs and cutoffs (a whole scale of
        // 1 or more here; `RpsRun::search_volume` refuses one whose products leave `Int4` and one
        // that changes Lambda).
        let scale = aux.scale;
        if !(scale >= 1.0 && scale.fract() == 0.0) {
            bail!(
                "a profile database of the scale {scale} (not a whole number of 1 or more) is not supported by LOSAT's {}",
                self.program
            );
        }
        // (The concatenated length is refused, or not, once for the run by its unit:
        // `RpsUnits::check_engine_lengths`.)
        // The consensus letters at the profiles' offsets (the traceback's query letters, read
        // as stored: `Blast_TracebackGetEncoding` gives `eBlastEncodingProtein`).
        Ok(())
    }

    /// The consensus letters of the volume's profiles at their offsets (`Consensus`).
    fn consensus(&self, volume: &Volume, files: &RpsVolumeFiles) -> Result<Vec<u8>> {
        let num_seqs = volume.num_oids();
        let concat = files.rows().len();
        let offsets = files.offsets();
        let mut consensus = vec![0u8; concat];
        if num_seqs > 0 {
            let block = volume.read_block(0, num_seqs, 0, 0).map_err(|error| {
                crate::algorithm::common::database::database_error(error, self.program)
            })?;
            for oid in 0..num_seqs as usize {
                let letters = block.protein(oid);
                let start = offsets[oid] as usize;
                consensus[start..start + letters.len()].copy_from_slice(letters);
            }
        }
        Ok(consensus)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:177-190
// ```c
//    case eRPSLookupTable:
//        {
//            BlastRPSLookupTable *lookup;
//            Int4 alphabet_size;
//            RPSLookupTableNew(rps_info, (BlastRPSLookupTable* *)
//                                         (&lookup_wrap->lut));
//
//            /* if the alphabet size from the RPS database is too
//               small, mask all unsupported query letters */
//            lookup = (BlastRPSLookupTable*)(lookup_wrap->lut);
//            alphabet_size = lookup->alphabet_size;
//            if (alphabet_size < BLASTAA_SIZE)
//                Blast_MaskUnsupportedAA(query, alphabet_size);
//            break;
//        }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1336-1346
// ```c
// Blast_MaskUnsupportedAA(BLAST_SequenceBlk* seq, Uint1 min_invalid)
// {
//     ...
//     for (i = 0; i < length; i++) {
//         if (sequence[i] >= min_invalid) {
//             sequence[i] = kProtMask;
//         }
//     }
// }
// ```
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
/// ```c++
///     while(input_it != input_list.end())
///     {
///         double final_evalue;
///         double input_evalue;
///
///         (*final_it)->GetNamedScore(CSeq_align::eScore_EValue, final_evalue);
///         (*input_it)->GetNamedScore(CSeq_align::eScore_EValue, input_evalue);
///
///         if(input_evalue == final_evalue)
///         {
///             //Pulling a trick here to keep the program flow simple
///             //Replace the final evalue with input bitscore and vice versa
///             (*final_it)->GetNamedScore(CSeq_align::eScore_BitScore, input_evalue);
///             (*input_it)->GetNamedScore(CSeq_align::eScore_BitScore, final_evalue);
///         }
///
///         if(input_evalue <  final_evalue)
///         {
///             ... (the whole group of the input's subject)
///             final_list.insert(final_it, start_input_it, input_it);
///         }
///         else
///         {
///             ... (final_it past the whole group of its subject)
///             if(final_it == final_list.end())
///             {
///                 final_list.insert(final_it, input_it, input_list.end());
///                 break;
///             }
///         }
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:155-160
/// ```c++
///             if(thread_results[d]->HasAlignments())
///             {
///                 CConstRef<CSeq_align_set>  thread_align_set = thread_results[d]->GetSeqAlign();
///                 if(align_set->IsEmpty())
///                 {
///                     align_set->Set().insert(align_set->Set().begin(),
///                                         thread_align_set->Get().begin(),
///                                         thread_align_set->Get().end());
/// ```
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
/// ```c++
/// #define SMALLEST_EVALUE 1.0e-180
/// ...
///     double evalue = (hsp->evalue < SMALLEST_EVALUE) ? 0.0 : hsp->evalue;
/// ```
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
mod tests {
    use super::{RpsUnit, RpsUnits};

    /// A unit of `total` letters (before the fallback) whose setup reads `setup` (`T`, `n`).
    fn unit(total: i64, setup: (i64, i32), request_count: i32, engine_size: i128) -> RpsUnit {
        RpsUnit {
            total,
            setup_total: setup.0,
            setup_count: setup.1,
            request_count,
            min_length: 10,
            engine_size,
            sorts_results: total > 0,
        }
    }

    fn units(units: Vec<RpsUnit>, adjusted: Option<(i64, i32)>, dbsize: i64) -> RpsUnits {
        let (db_length, db_num_seqs) = adjusted.unwrap_or((dbsize, 0));
        let order = (0..units.len()).collect();
        RpsUnits {
            units,
            adjusted,
            db_length,
            db_num_seqs,
            order,
        }
    }

    fn blosum62_gapped() -> crate::stats::KarlinParams {
        crate::stats::KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        }
    }

    /// Model M2 of line RG D1a-4: `D` and `M` above 0 replace each unit's length and count, at
    /// 0 or below each unit reads its own (and, at 0, its own Gumbel length); `-dbsize` replaces
    /// the lengths, never the counts.
    #[test]
    fn each_unit_reads_its_own_numbers_under_the_options() {
        let volume = unit(190, (190, 1), 1, 191);
        for (adjusted, eff, gumbel) in [
            (Some((560, 3)), (560, 3), 560),
            (Some((5000, 7)), (5000, 7), 5000),
            (Some((0, 0)), (190, 1), 190),
            (Some((0, -5)), (190, 1), 190),
            (
                Some((-9_223_372_036_854_775_807, 3)),
                (190, 3),
                -9_223_372_036_854_775_807,
            ),
        ] {
            let all = units(vec![volume; 3], adjusted, 0);
            assert_eq!(all.eff_lengths(&volume), eff, "{adjusted:?}");
            assert_eq!(all.gumbel_length(&volume), gumbel, "{adjusted:?}");
        }
        let dbsize = units(vec![volume; 3], None, 100_000);
        assert_eq!(dbsize.eff_lengths(&volume), (100_000, 1));
        assert_eq!(dbsize.gumbel_length(&volume), 100_000);
        let negative = units(vec![volume; 3], None, -5);
        assert_eq!(negative.eff_lengths(&volume), (190, 1));
        assert_eq!(negative.gumbel_length(&volume), -5);
        // One volume, LENGTH 0 and NSEQ 600: the setup's fallback to OID 0 and one sequence.
        let fallback = unit(0, (190, 1), 600, 600);
        let one = units(vec![fallback], None, 0);
        assert_eq!(one.eff_lengths(&fallback), (190, 1));
        assert_eq!(one.gumbel_length(&fallback), 0);
    }

    /// Model M6: the request's search space is the last unit's, with its total and count before
    /// the fallback: 0 for a total of 0, the query length for a negative one.
    #[test]
    fn the_request_space_is_the_last_units_before_the_fallback() {
        let gapped = blosum62_gapped();
        let first = unit(190, (190, 1), 1, 191);
        let last = unit(180, (180, 1), 1, 181);
        let two = units(vec![first, last], Some((0, 0)), 0);
        assert_eq!(
            two.request_space(208, &gapped, "RPSBLAST").unwrap(),
            two.search_space(&last, 208, &gapped, 0, "RPSBLAST")
                .unwrap()
        );
        assert_ne!(
            two.request_space(208, &gapped, "RPSBLAST").unwrap(),
            two.search_space(&first, 208, &gapped, 0, "RPSBLAST")
                .unwrap()
        );
        let zero = units(vec![unit(0, (190, 1), 600, 600)], None, 0);
        assert_eq!(zero.request_space(208, &gapped, "RPSBLAST").unwrap(), 0);
        let negative = units(
            vec![unit(-9_223_372_036_854_775_807, (190, 1), 600, 601)],
            None,
            0,
        );
        assert_eq!(
            negative.request_space(208, &gapped, "RPSBLAST").unwrap(),
            208
        );
        assert_eq!(negative.request_space(0, &gapped, "RPSBLAST").unwrap(), 0);
    }

    // prelim_stage.cpp:219-221 (quoted at `RpsUnits::prelim_lengths`): the preliminary stage reads
    // SplitQuery's space, whose lengths are the setup's except for a negative total (review A of
    // D1a-4, MF1): `LENGTH` 2^63 + 1 and `NSEQ` 600 over one volume, SplitQuery's 208 (1 times the
    // query's 208 letters less no adjustment) where the setup reads OID 0's 190 letters and 1.
    #[test]
    fn the_preliminary_lengths_differ_from_the_setups_for_a_negative_total_only() {
        let negative = unit(-9_223_372_036_854_775_807, (190, 1), 600, 601);
        let one = |unit: RpsUnit| units(vec![unit], None, 0);
        assert_eq!(
            one(negative).prelim_lengths(&negative, 0),
            Some((-9_223_372_036_854_775_807, 600))
        );
        // `-searchsp`: both stages read it.
        assert_eq!(one(negative).prelim_lengths(&negative, 1_000_000), None);
        // `LENGTH 0`: SplitQuery sets no space, the setup computes its own.
        let zero = unit(0, (190, 1), 600, 600);
        assert_eq!(one(zero).prelim_lengths(&zero, 0), None);
        // A positive total: the same lengths (also `STATS_NSEQ`'s count, and under D and M).
        let plain = unit(560, (560, 7), 7, 563);
        assert_eq!(one(plain).prelim_lengths(&plain, 0), None);
        let volume = unit(190, (190, 1), 1, 191);
        let several = units(vec![volume, volume], Some((-5, -3)), 0);
        assert_eq!(several.prelim_lengths(&volume, 0), None);
        let several = units(vec![volume, volume], Some((5000, 7)), 0);
        assert_eq!(several.prelim_lengths(&volume, 0), None);
    }

    // blast_setup.c:826-846 (quoted at `stats::setup::calc_eff_lengths`): with `-dbsize` each
    // volume's space reads the given length and its own sequence count; an effective database
    // length of 0 or below is 1; a product past `Int8` is refused.
    #[test]
    fn a_units_search_space_reads_its_own_sequence_count() {
        use crate::stats::length_adjustment::compute_length_adjustment_ncbi;
        let gapped = blosum62_gapped();
        for (query, db_length, seqs) in
            [(302i64, 5_000i64, 1i32), (302, 5_000, 3), (40, 100_000, 7)]
        {
            let volume = unit(190, (190, seqs), seqs, 191);
            let adjustment =
                compute_length_adjustment_ncbi(query, db_length, i64::from(seqs), &gapped)
                    .length_adjustment;
            assert_eq!(
                units(vec![volume; 3], None, db_length)
                    .search_space(&volume, query, &gapped, 0, "RPSBLAST")
                    .unwrap(),
                (db_length - i64::from(seqs) * adjustment) * (query - adjustment),
                "{query} {db_length} {seqs}"
            );
        }
        let volume = unit(190, (190, 1), 1, 191);
        assert_eq!(
            units(vec![volume], None, 0)
                .search_space(&volume, 302, &gapped, 1_000_000, "RPSBLAST")
                .unwrap(),
            1_000_000
        );
        assert!(units(vec![volume; 3], None, i64::MAX)
            .search_space(&volume, 302, &gapped, 0, "RPSBLAST")
            .unwrap_err()
            .to_string()
            .contains("overflows NCBI BLAST+'s 64-bit effective search space"));
    }

    /// The strategy's `DbLength`: `D` of any sign where `x_AdjustDbSize` runs; -2^63 refused.
    #[test]
    fn the_strategy_db_length_is_d_and_minus_two_to_the_63_is_refused() {
        let volume = unit(190, (190, 1), 1, 191);
        for length in [0, -1, 5000, -9_223_372_036_854_775_807] {
            assert_eq!(
                units(vec![volume; 3], Some((length, 3)), 0)
                    .strategy_db_length("RPSBLAST")
                    .unwrap(),
                Some(length)
            );
        }
        assert_eq!(
            units(vec![volume], None, 0)
                .strategy_db_length("RPSBLAST")
                .unwrap(),
            None
        );
        let refusal = units(vec![volume; 3], Some((i64::MIN, 3)), 0)
            .strategy_db_length("RPSTBLASTN")
            .unwrap_err()
            .to_string();
        assert!(
            refusal.contains("in the search strategy's DbLength"),
            "{refusal}"
        );
        assert!(refusal.ends_with("LOSAT's RPSTBLASTN"), "{refusal}");
    }

    /// The setup's division by 0 without the Gumbel block, and the concatenated lengths.
    #[test]
    fn counts_of_zero_and_concatenated_lengths_are_refused() {
        let zero = units(vec![unit(560, (560, 0), 0, 560)], None, 0);
        assert!(zero.check_setup_counts(true, "RPSBLAST").is_ok());
        assert!(zero.check_setup_counts(false, "RPSBLAST").is_err());
        let check = |engine_size: i128| {
            units(vec![unit(560, (560, 3), 3, engine_size)], None, 0)
                .check_engine_lengths(40, "RPSBLAST")
                .map_err(|error| error.to_string())
        };
        assert!(check(563).is_ok());
        assert!(check(1003).is_ok());
        // LENGTH 0 and NSEQ -300 over the 563 rows of three profiles: searched (an `Lc` below the
        // rows, `RpsRun::concat_length`), also a negative one (LENGTH 0, NSEQ -5).
        assert!(check(3).is_ok());
        assert!(check(260).is_ok());
        assert!(check(-2).is_ok());
        // LENGTH 2^63 + 1 and NSEQ 600: `(Int4)` keeps 601.
        assert!(check(-9_223_372_036_854_775_807 + 600).is_ok());
        // The diagonal table: NCBI's loop ends at 2^30 (`Lc` plus the window 40 of exactly 2^30,
        // also with `(Int4)` of a negative sum) and never past it (LENGTH 2^30).
        assert!(check((1 << 30) - 40).is_ok());
        assert!(check((1 << 30) - 40 - (1i128 << 32)).is_ok());
        assert!(check((1 << 30) - 39)
            .unwrap_err()
            .contains("never stops growing"));
        assert!(check(1_073_741_827)
            .unwrap_err()
            .contains("of 1073741827 concatenated positions"));
        // LENGTH 2^31 - 4 and NSEQ 3: `Lc` is INT4_MAX, whose sum with the window overflows NCBI's
        // `Int4`; the widened sum is above 2^30.
        assert!(check(2_147_483_647)
            .unwrap_err()
            .contains("never stops growing"));
        // LENGTH 2^31: above INT4_MAX (the double free comes first in NCBI's engine).
        assert!(check(2_147_483_651)
            .unwrap_err()
            .contains("above 2147483647"));
        // The window of the one-hit finder (`-window_size 0`): the sum is `Lc` alone.
        let one_hit = |engine_size: i128| {
            units(vec![unit(560, (560, 3), 3, engine_size)], None, 0)
                .check_engine_lengths(0, "RPSBLAST")
                .is_ok()
        };
        assert!(one_hit(1 << 30));
        assert!(!one_hit((1 << 30) + 1));
    }
}
