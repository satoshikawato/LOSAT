//! The volumes of a profile database as `rpsblast` and `rpstblastn` search them (series PR,
//! items PS11b and PS12a): NCBI's order of the volumes and their first ordinal ids, the set-up of
//! one volume (its files, LOSAT's refusals of what NCBI's search does not reproducibly run, the
//! consensus letters of its profiles), the masking of the letters that a profile has no column
//! for, and the merge of the volumes' results. One owner for both programs (review A of PS12a, N3).

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
}

/// Whether the volumes of an RPS search have their own effective search spaces: `-dbsize`
/// (`dbsize` above 0) with several volumes and no `-searchsp` above 0 (`x_AdjustDbSize` returns
/// early, so each volume's statistics read the given length and the volume's own sequence count).
pub(crate) fn per_volume_search_spaces(dbsize: i64, searchsp: i64, volumes: usize) -> bool {
    dbsize > 0 && searchsp == 0 && volumes >= 2
}

/// One volume set up for a search: its files and the consensus letters of its profiles at their
/// offsets (the traceback's letters and the reports' identities).
pub(crate) struct RpsVolume {
    pub(crate) files: RpsVolumeFiles,
    pub(crate) consensus: Vec<u8>,
}

/// The set-up of the volumes of one run: the program of LOSAT's refusals, the word window (the
/// diagonal table's bound) and the matrix and gap costs of the first volume set up.
pub(crate) struct RpsVolumeSetup {
    program: &'static str,
    window: i32,
    first_costs: Option<(String, i32, i32)>,
}

impl RpsVolumeSetup {
    pub(crate) fn new(program: &'static str, window: i32) -> Self {
        Self {
            program,
            window,
            first_costs: None,
        }
    }

    /// `volume` set up for a search, with its frequency ratios when `cbs` (composition-based
    /// statistics), after LOSAT's checks (`check`): `Ok(Err(error))` is NCBI's own error of the
    /// files, which the caller reports where NCBI's search raises it; `Err` is LOSAT's refusal.
    pub(crate) fn open(
        &mut self,
        volume: &Volume,
        cbs: bool,
    ) -> Result<std::result::Result<RpsVolume, RpsVolumeError>> {
        let files = match self.files(volume, cbs)? {
            Ok(files) => files,
            Err(error) => return Ok(Err(error)),
        };
        self.check_files(volume, &files)?;
        let consensus = self.consensus(volume, &files)?;
        Ok(Ok(RpsVolume { files, consensus }))
    }

    /// LOSAT's checks of `volume` alone, once per run before the first byte of output
    /// (rpstblastn, which maps the volumes again for each query batch: `reopen`).
    pub(crate) fn check(
        &mut self,
        volume: &Volume,
        cbs: bool,
    ) -> Result<std::result::Result<(), RpsVolumeError>> {
        let files = match self.files(volume, cbs)? {
            Ok(files) => files,
            Err(error) => return Ok(Err(error)),
        };
        self.check_files(volume, &files)?;
        Ok(Ok(()))
    }

    /// `volume`, which `check` has checked in this run, set up again for a search.
    pub(crate) fn reopen(
        &self,
        volume: &Volume,
        cbs: bool,
    ) -> Result<std::result::Result<RpsVolume, RpsVolumeError>> {
        let files = match self.files(volume, cbs)? {
            Ok(files) => files,
            Err(error) => return Ok(Err(error)),
        };
        let consensus = self.consensus(volume, &files)?;
        Ok(Ok(RpsVolume { files, consensus }))
    }

    /// The volume's files (and frequency ratios when `cbs`).
    fn files(
        &self,
        volume: &Volume,
        cbs: bool,
    ) -> Result<std::result::Result<RpsVolumeFiles, RpsVolumeError>> {
        let num_seqs = volume.num_oids();
        let lengths = (0..num_seqs)
            .map(|oid| volume.protein_length(oid))
            .collect::<std::result::Result<Vec<u32>, _>>()
            .map_err(|error| {
                crate::algorithm::common::database::database_error(error, self.program)
            })?;
        let mut files = match RpsVolumeFiles::open_for(volume.base(), &lengths, self.program) {
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

    /// The consensus letters of the volume's profiles at their offsets.
    fn consensus(&self, volume: &Volume, files: &RpsVolumeFiles) -> Result<Vec<u8>> {
        let num_seqs = volume.num_oids();
        let concat = files.rows().len();
        let offsets = files.offsets();
        let mut consensus = vec![0u8; concat];
        if num_seqs > 0 {
            let block = volume.read_block(0, num_seqs, 0).map_err(|error| {
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
