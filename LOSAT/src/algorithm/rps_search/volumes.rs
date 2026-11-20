//! The volumes of a profile database as `rpsblast` and `rpstblastn` search them (series PR,
//! items PS11b and PS12a): NCBI's order of the volumes and their first ordinal ids, the set-up of
//! one volume (its files, LOSAT's refusals of what NCBI's search does not reproducibly run, the
//! consensus letters of its profiles), the masking of the letters that a profile has no column
//! for, and the merge of the volumes' results. One owner for both programs (review A of PS12a, N3).

use anyhow::{bail, Result};

use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::seqdb::rps::{canonical_volume_path, RpsOpenError, RpsVolumeFault, RpsVolumeFiles};

/// NCBI's own errors of a volume's files, which its search of a query batch raises (`BLAST
/// engine error`, exit status 3).
pub(crate) enum RpsVolumeError {
    /// NCBI's text of a missing, short or unreadable file.
    Files(String),
    /// An `.aux` without a matrix name (`RpsVolumeFault::NoMatrix`): NCBI's score block fails for
    /// each query of the batch ("Unknown error code -1").
    NoMatrix,
}
use crate::seqdb::{DbStatistics, SeqDb, Volume};

/// LOSAT's refusal of a profile database whose total length (`BlastSeqSrcGetTotLen`, NCBI's `Int8`)
/// is 0 or negative: an alias file's `LENGTH 0` or a `LENGTH` of 2^63 or more (no list reaches it:
/// NCBI's RPS programs take no list option, blast_args.cpp:2241-2251, and the alias keys of lists
/// are refused before). NCBI's RPS search has no loop over subjects (blast_engine.c:1384-1393) and
/// opens each volume alone, whose statistics then read the volume's own length and count (probe of
/// D1a-3 fix 6: a footer's "Effective search space used" and the next request's
/// `EffectiveSearchSpace` differ from a plain reading of the negative total); LOSAT has not traced
/// these totals (item D1a-4). One owner of the refusal and its text for `rpsblast` and `rpstblastn`
/// (`program`), so neither program reaches `common::database` for it.
pub(crate) fn untraced_total_length(statistics: &DbStatistics, program: &str) -> Result<()> {
    if statistics.engine_db_length <= 0 {
        return Err(crate::algorithm::common::database::refusal(
            "an alias file's LENGTH of 0 or of 2^63 or more (a total length of 0 or less as NCBI's Int8) over an RPS-BLAST database",
            program,
        ));
    }
    Ok(())
}

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
        // The first ordinal id of each volume (the database's order), and the volumes in the order
        // of their paths as strings (a name given twice is searched once).
        let mut starts = Vec::with_capacity(db_volumes.len());
        let mut next_start = 0u32;
        for volume in db_volumes {
            starts.push(next_start);
            next_start += volume.num_oids();
        }
        let mut order: Vec<usize> = (0..db_volumes.len()).collect();
        order.sort_by(|&a, &b| {
            crate::seqdb::SeqDb::compare_volume_paths(db_volumes[a].base(), db_volumes[b].base())
        });
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
/// remote option), then each volume's set-up the matrix and the gap costs of its `.aux`
/// (`CreateRpsStructures`; the volumes' are the same, `RpsVolumeSetup`), also in a batch without a
/// valid query. `costs` are the volumes' matrix and gap costs as their files' set-up read them
/// (`RpsVolumeSetup::costs`), `None` when the files do not open: the batch's search fails there with
/// NCBI's error and no request follows. NCBI's threads over the volumes set them up on clones
/// (`RunThreadedSearch`), whose search spaces never reach the handle; LOSAT writes the requests of
/// NCBI's `-num_threads 1` at every thread count (DECISIONS F1b-2 of line RG).
pub(crate) fn search_sets(
    db: &SeqDb,
    costs: Option<(String, i32, i32)>,
    searchsp: Option<i64>,
    dbsize: Option<i64>,
) -> Vec<crate::report::archive::options::Param> {
    use crate::report::archive::options::{Param, ParamValue};
    let volumes = RpsVolumeSet::of(db);
    let mut sets = Vec::new();
    if volumes.order.len() > 1 && searchsp.unwrap_or(0) == 0 && dbsize.unwrap_or(0) == 0 {
        // `check_alias_keys` refuses `STATS_TOTLEN`, so the length is `GetTotalLength`, the `Uint8`
        // that `SetDbLength(Int8)` takes as NCBI's `Int8`: `DbStatistics::engine_db_length`
        // (positive here: `untraced_total_length` refuses 0 and the negative totals of an alias
        // `LENGTH` of 2^63 or more).
        sets.push(Param {
            name: "DbLength",
            value: ParamValue::BigInteger(db.statistics().engine_db_length),
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

/// Whether each volume of an RPS search reads its own statistics: with several volumes, a
/// `-searchsp` or a `-dbsize` makes `x_AdjustDbSize` return before it sets the database's length
/// and count, so that each volume's search reads its own sequence count and, without `-dbsize`,
/// its own total length (`real_db_length`); the reports print volume 0's (`first`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-727
/// ```c
///    if (eff_len_options->db_length > 0)
///       db_length = eff_len_options->db_length;
///    else
///       db_length = eff_len_params->real_db_length;
///    ...
///    if (eff_len_options->dbseq_num > 0)
///       db_num_seqs = eff_len_options->dbseq_num;
///    else
///       db_num_seqs = eff_len_params->real_num_seqs;
/// ```
pub(crate) fn per_volume_lengths(dbsize: i64, searchsp: i64, volumes: usize) -> bool {
    volumes >= 2 && (dbsize > 0 || searchsp != 0)
}

/// The total length that a volume's statistics read without `-dbsize` (`real_db_length` of the
/// volume's own database): its letters.
pub(crate) fn volume_db_length(volume: &Volume) -> i64 {
    i64::try_from(volume.index().volume_length).unwrap_or(i64::MAX)
}

/// Whether the volumes of an RPS search have their own effective search spaces: `-dbsize`
/// (`dbsize` above 0) with several volumes and no `-searchsp` above 0 (`x_AdjustDbSize` returns
/// early, so each volume's statistics read the given length and the volume's own sequence count).
pub(crate) fn per_volume_search_spaces(dbsize: i64, searchsp: i64, volumes: usize) -> bool {
    dbsize > 0 && searchsp == 0 && volumes >= 2
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

/// The set-up of the volumes of one run: the program of LOSAT's refusals, the word window (the
/// diagonal table's bound) and the matrix and gap costs of the first volume set up.
#[derive(Debug, PartialEq)]
pub(crate) struct RpsVolumeSetup {
    program: &'static str,
    window: i32,
    first_costs: Option<(String, i32, i32)>,
    /// The volumes (`canonical_volume_path`) whose cells `check` or `open` found not in groups of
    /// buckets in order and regrouped (`RpsVolumeFiles::is_regrouped`; series PR, item PS12m):
    /// `reopen` regroups them again without reading every hit. This set-up is the one owner of the
    /// fact: a run that reopens volumes keeps the set-up that checked them.
    regrouped: std::collections::BTreeSet<String>,
    /// The `-db` value of a database of one volume, which NCBI's error of a missing `.freq`
    /// names (`with_database`); `None`: each volume's path.
    database_name: Option<String>,
}

impl RpsVolumeSetup {
    pub(crate) fn new(program: &'static str, window: i32) -> Self {
        Self {
            program,
            window,
            first_costs: None,
            regrouped: std::collections::BTreeSet::new(),
            database_name: None,
        }
    }

    /// The set-up of the volumes of the database `name` (the `-db` value as NCBI keeps it) of
    /// `num_volumes` volumes: NCBI searches a database of one volume by its name, so that the
    /// error of a missing `.freq` names it as given, and each volume of several by its path.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:341-346, 380-386
    /// ```c++
    /// 	CSeqDB::FindVolumePaths(db, CSeqDB::eProtein, m_rps_databases, NULL, true, true);
    /// 	m_num_of_dbs = static_cast<unsigned int>(m_rps_databases.size());
    /// 	...
    /// 		if(1 == m_num_of_dbs)
    /// 		{
    /// 			return s_RunLocalRpsSearch(m_db_name, *m_query_vector, m_opt_handle);
    /// ```
    pub(crate) fn with_database(mut self, name: &str, num_volumes: usize) -> Self {
        self.database_name = (num_volumes == 1).then(|| name.to_string());
        self
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
        let concat = files.rows().len();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1129-1133
        // ```c
        //     num_db_seqs = BlastSeqSrcGetNumSeqs(seq_src);
        //     dbsize = BlastSeqSrcGetTotLen(seq_src) + num_db_seqs;
        //     if (dbsize > INT4_MAX)
        //       return -3;
        // ```
        if concat as u64 + self.window as u64 > (i32::MAX / 2) as u64 {
            bail!(
                "a profile volume of {concat} positions is not supported by LOSAT's {}",
                self.program
            );
        }
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
    use super::untraced_total_length;
    use crate::seqdb::DbStatistics;

    fn with_total(engine_db_length: i64) -> DbStatistics {
        DbStatistics {
            setup_total: engine_db_length,
            setup_count: 3,
            gumbel_db_length: engine_db_length,
            engine_db_length,
            request_count: 3,
            min_length: 10,
            average: 100,
        }
    }

    /// The RPS search refuses a total length of 0 (an alias `LENGTH 0`) and the negative totals of
    /// an alias `LENGTH` of 2^63 or more (`Int8` of `GetTotalLength`), and searches a positive one.
    #[test]
    fn an_rps_search_refuses_totals_of_zero_or_less() {
        for total in [0, i64::MIN, i64::MIN + 1, -1] {
            let refusal = untraced_total_length(&with_total(total), "rpsblast").unwrap_err();
            assert_eq!(
                refusal.to_string(),
                "an alias file's LENGTH of 0 or of 2^63 or more (a total length of 0 or less as \
                 NCBI's Int8) over an RPS-BLAST database; this database is not supported by \
                 LOSAT's rpsblast"
            );
        }
        assert!(untraced_total_length(&with_total(1), "rpsblast").is_ok());
        assert!(untraced_total_length(&with_total(i64::MAX), "rpstblastn").is_ok());
    }
}
