//! The search of a run's queries against one profile volume: every (query, volume) pair alone,
//! each valid context of each query (a protein query's one, a translated query's frames) a unit
//! of the run's pool on a worker's scratch, the query's traceback run by the worker that searches
//! its last context, the results in query order (series PR, item PS12s, D2). The outfmt 0 K of
//! each query batch is computed from the queries' marks after (`batch_gapped_k`).
//!
//! NCBI searches the queries of a batch one after another on one diagonal table, then traces the
//! lists of the whole batch; a query's results depend on the volume, the options and the query
//! only (the table's entries of earlier queries and contexts read as empty, `wordfinder.rs`), so
//! LOSAT's bytes are the same at every thread count.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1161-1188

use anyhow::Result;

use super::prelim::{RpsPrelim, RpsPrelimScratch};
use super::redo::RpsRedoParams;
use super::scan::{RpsBlockPool, RpsTable};
use super::traceback::{
    filter_results, formatter_identities, RpsTraceScratch, RpsTraceback, RpsTraced,
};
use super::wordfinder::{RpsWordParams, RpsWordScratch};
use super::{RpsPastRows, RpsRows, RPS_K_MULT};
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::pssm_search::evalue::EvalueStats;
use crate::algorithm::pssm_search::redo1::{RedoRun, SCALING_FACTOR};
use crate::seqdb::rps::RpsVolumeFiles;
use crate::stats::{BlastGumbelBlk, KarlinParams};

/// The options of an RPS-BLAST search that blastp's driver does not apply itself: `-max_hsps`
/// and `-qcov_hsp_perc` (the RPS traceback filters with the length of the first query of a
/// batch, `traceback::filter_results`), and the composition-based statistics mode that the
/// command line asks for (the mode the search runs is blastp's `comp_based_stats`, which
/// `rps_hooks::check_freq_ratio_file` sets to 0 where the database has no frequency ratios).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RpsSearch {
    pub(crate) max_hsps_per_subject: usize,
    pub(crate) query_cov_hsp_perc: f64,
    pub(crate) requested_cbs: u8,
    /// `-num_threads 0`: one thread per volume of the database (NCBI's "auto").
    pub(crate) auto_threads: bool,
    /// `-xdrop_ungap`, `-xdrop_gap` and `-xdrop_gap_final` in bits as given; `None`: NCBI's RPS
    /// defaults (7, 15 and 25 bits, which `rps_hooks::rps_search` reads).
    pub(crate) x_drop_ungapped_bits: Option<f64>,
    pub(crate) x_drop_gapped_bits: Option<f64>,
    pub(crate) x_drop_final_bits: Option<f64>,
    /// `-sum_stats true`, refused where the search starts (`rps_hooks::rps_search`).
    pub(crate) sum_stats: bool,
    /// `-use_sw_tback`: the redo of composition-based statistics 1 to 3 is NCBI's
    /// Smith-Waterman redo (`RpsRedoParams::smith_waterman`), of rpstblastn's frames too (series
    /// PR, item PS12c); without them NCBI's traceback does not read it.
    pub(crate) use_sw_tback: bool,
    /// Whether this RPS search is DELTA-BLAST's domain search (`CDeltaBlast::x_FindDomainHits`):
    /// refusals name `"DELTABLAST"`, the one-time `.loo` check releases its mapped pages before
    /// the search so single-query iterations do not hold unused pages, and concurrent diagonal
    /// tables are bounded so peak RSS stays at or below NCBI's single-threaded domain search.
    pub(crate) domain_search: bool,
}

/// The contexts of a translated query.
pub(crate) const NUM_FRAMES: usize = 6;

impl RpsSearch {
    /// The program name shown in LOSAT's refusals for this RPS search.
    pub(crate) fn program(&self) -> &'static str {
        if self.domain_search {
            "DELTABLAST"
        } else {
            "RPSBLAST"
        }
    }
    /// The X-drops in bits, `-xdrop_ungap`, `-xdrop_gap` and `-xdrop_gap_final` as given or the
    /// RPS options handle's defaults (7, 15 and 25 bits), and the gap trigger (22 bits): the
    /// `RpsRun` fields of both RPS programs.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:76-89
    /// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:122,128,137,144
    pub(crate) fn x_drop_bits(&self) -> RpsXDropBits {
        RpsXDropBits {
            ungapped: self.x_drop_ungapped_bits.unwrap_or(7.0),
            gapped: self.x_drop_gapped_bits.unwrap_or(15.0),
            final_: self.x_drop_final_bits.unwrap_or(25.0),
            gap_trigger: 22.0,
        }
    }
}

/// `RpsSearch::x_drop_bits`.
pub(crate) struct RpsXDropBits {
    pub(crate) ungapped: f64,
    pub(crate) gapped: f64,
    pub(crate) final_: f64,
    pub(crate) gap_trigger: f64,
}

/// One user query as the search reads it: its contexts, one for a protein query (rpsblast) and
/// the six frames of a nucleotide query (rpstblastn), which the reversed search scans as its
/// subjects.
pub(crate) struct RpsQuery<'a> {
    /// The contexts in NCBI's order (frame 0; or +1, +2, +3, -1, -2, -3).
    pub(crate) frames: Vec<RpsFrame<'a>>,
    /// The ungapped block of the first valid query of the query's batch: `kbp_psi[0]` of the
    /// score block that the composition-based redo copies, which stays where the profile has no
    /// score frequencies.
    pub(crate) batch_kbp: KarlinParams,
}

/// One context of a user query: a protein query's letters (frame 0) or one frame of a
/// nucleotide query (NCBI's `subject->frame`, `subject->sequence` and `subject->length` of the
/// frame loop).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-813
#[derive(Debug, Clone, Copy)]
pub(crate) struct RpsFrame<'a> {
    /// 0 for a protein query; +1, +2, +3, -1, -2 or -3.
    pub(crate) frame: i8,
    /// The letters between two sentinels as searched (SEG, lower-case masks,
    /// `Blast_MaskUnsupportedAA`); a frame without letters is two sentinels.
    pub(crate) aa_seq: &'a [u8],
    /// The letters as read (between two sentinels) where the searched ones differ: the reports'
    /// identities count on them (`traceback::formatter_identities`).
    pub(crate) aa_seq_nomask: Option<&'a [u8]>,
    /// Whether its Karlin-Altschul parameters exist (an invalid context gets no HSP: its cutoffs
    /// are `INT4_MAX`, blast_parameters.c:916-919).
    pub(crate) valid: bool,
    /// `kbp_std[context]`: the context's composition with BLOSUM62.
    pub(crate) kbp_std: KarlinParams,
    /// `query_info->contexts[context].eff_searchsp` in this volume (NCBI's `Int8`): what the
    /// e-values and the hit cutoff read without the Gumbel block (`EvalueStats`), and the
    /// Smith-Waterman redo with it (series PR, item PS11c-3). 0 for an invalid context.
    pub(crate) search_space: i64,
    /// The length of the first query of the query's batch at this context (`-qcov_hsp_perc`).
    pub(crate) first_query_length: usize,
}

impl RpsFrame<'_> {
    /// The letters without the sentinels (`subject->length`).
    pub(crate) fn length(&self) -> usize {
        self.aa_seq.len().saturating_sub(2)
    }

    /// What the context's e-values read (`EvalueStats`): the Gumbel block, if any, and the
    /// context's search space.
    pub(super) fn stats<'g>(&self, gumbel: Option<&'g BlastGumbelBlk>) -> EvalueStats<'g> {
        EvalueStats::new(gumbel, self.search_space)
    }
}

impl RpsQuery<'_> {
    /// The context of `frame` (`BLAST_FrameToContext`; 0 for frame 0).
    pub(crate) fn context(&self, frame: i32) -> &RpsFrame<'_> {
        &self.frames[crate::core::blast_util::frame_to_context(frame)]
    }

    /// The length of the query's translation buffer (NCBI's `one_query->length`; a protein
    /// query's length): the last context's offset plus its length, where each context starts
    /// one sentinel after the previous one, an empty context taking no byte.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:72-77
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:379-381
    pub(crate) fn translation_length(&self) -> usize {
        let (last, earlier) = self.frames.split_last().expect("a query has a context");
        earlier
            .iter()
            .map(|frame| match frame.length() {
                0 => 0,
                length => length + 1,
            })
            .sum::<usize>()
            + last.length()
    }

    /// `stat_length` of the preliminary e-values: the length of the last context that has
    /// letters (a protein query's length; frame -3 with both strands, +3 with `-strand plus`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:813 (cited at
    /// `RpsFrame`)
    pub(crate) fn stat_length(&self) -> usize {
        self.frames
            .iter()
            .rev()
            .map(RpsFrame::length)
            .find(|&length| length > 0)
            .unwrap_or(0)
    }
}

/// What the search of every query of a volume reads besides the volume.
pub(crate) struct RpsRun<'a> {
    pub(crate) queries: Vec<RpsQuery<'a>>,
    /// `kbp_gap_std`: BLOSUM62 11/1's table values.
    pub(crate) gapped: KarlinParams,
    /// `sbp->gbp`: none with `-searchsp` (at any value) or `OLD_FSC`.
    pub(crate) gumbel: Option<&'a BlastGumbelBlk>,
    /// `avg_subject_length` of the hit cutoff (`BLAST_SEQSRC_MINLENGTH`).
    pub(crate) cutoff_subject_length: i32,
    /// `-evalue`.
    pub(crate) evalue: f64,
    /// `hit_params->prelim_evalue` and the e-value of the hit cutoff: `-evalue`, five times it
    /// with composition-based statistics 2 and 3 and the Gumbel block (`cbs_stretch`,
    /// blast_parameters.c:937).
    pub(crate) prelim_evalue: f64,
    /// The composition-based statistics mode the search runs (0 to 3).
    pub(crate) cbs: u8,
    /// BLOSUM62's score block for the redo's matrix information (modes 1 to 3).
    pub(crate) sbp: Option<&'a crate::core::psi::PsiScoreBlk>,
    pub(crate) hitlist_size: usize,
    pub(crate) prelim_hitlist_size: usize,
    pub(crate) window: i32,
    /// The concatenated length `Lc` of the volume's search unit (`RpsUnit::concat_length`,
    /// `concat_db.length`), NCBI's "query" length of the search: the length of its diagonal
    /// table, the bound of the ungapped extensions to the right and the end of the gapped stage's
    /// interval tree. The volume's rows, except where an alias file's `LENGTH` and `NSEQ` over
    /// one volume move it (below the rows, also negative, or above them).
    pub(crate) concat_length: i32,
    /// Whether the traceback sorts each query's list again after its filters
    /// (`RpsUnit::sorts_results`: the unit's total length above 0; `finish_volume_list`).
    pub(crate) sorts_results: bool,
    /// The search spaces of the preliminary stage's hit cutoff (without the Gumbel block), one for
    /// each context of each query (query after query, `RpsQuery::frames` each), where they are not
    /// the contexts' `RpsFrame::search_space`: `SplitQuery_SetEffectiveSearchSpace`'s spaces over
    /// one volume whose total is negative (`RpsUnits::prelim_lengths`); `None` elsewhere.
    pub(crate) prelim_spaces: Option<&'a [i64]>,
    /// The X-drops in bits (`word_options->x_dropoff`, `ext_options->gap_x_dropoff` and
    /// `gap_x_dropoff_final`: `-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final` or NCBI's RPS
    /// defaults) and the gap trigger in bits (`kOptions->gap_trigger`).
    pub(crate) x_drop_ungapped_bits: f64,
    pub(crate) x_drop_gapped_bits: f64,
    pub(crate) x_drop_final_bits: f64,
    pub(crate) gap_trigger_bits: f64,
    pub(crate) options: &'a RpsSearch,
    /// The program of LOSAT's refusals (`RPSBLAST`, `RPSTBLASTN`).
    pub(crate) program: &'static str,
}

/// The refusal of a database whose scale makes a scaled score leave NCBI's `Int4` (gap costs,
/// X-drops, cutoffs): NCBI's products overflow there.
fn scale_overflow(scale: f64, program: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "a profile database of the scale {scale} (its scaled gap costs, X-drops or score cutoffs leave 32-bit integers) is not supported by LOSAT's {program}"
    )
}

/// The refusal of an X-drop in bits whose raw score (the bits over Lambda, converted to `Int4`,
/// times the database's scale where NCBI multiplies) leaves NCBI's `Int4`, or whose quotient by
/// the scaled gap extension makes NCBI's kernels' sizes leave it (`kernel_x_drop_fits`): NCBI's
/// conversion of the double, or its arithmetic, is undefined there (DECISIONS 17 N1/N2, line PS).
/// `given` is the option's value as the command line gives it; an X-drop of NCBI's defaults gets
/// the database's refusal (`scale_overflow`), as only an extreme scale makes it leave 32 bits.
fn x_drop_overflow(
    option: &str,
    given: Option<f64>,
    bits: f64,
    scale: f64,
    program: &str,
) -> anyhow::Error {
    if given.is_none() {
        return scale_overflow(scale, program);
    }
    anyhow::anyhow!(
        "an {option} of {bits:?} bits, whose raw X-drop at the profile database's scale {scale} leaves NCBI BLAST+'s 32-bit integers (its conversion or the sizes of its alignment kernels overflow: undefined behaviour), is not supported by LOSAT's {program}"
    )
}

// The limit of NCBI's gapped kernels' sizes (`kernel_x_drop_fits`, its `Int4` arithmetic on the X-drop
// over the scaled gap extension) is the protein programs' (`core::x_drop`). Reachable here only
// where the scaled gap extension is 1 or 2 (scale 1 or 2 with BLOSUM62 11/1): the quotient is at
// most `INT4_MAX` over `scale` otherwise.
use crate::core::x_drop::kernel_x_drop_fits;

/// A query's ungapped X-drop (`x_dropoff_init`) from the X-drop in bits `bits` (above 0,
/// `validate_protein_options`), the query's ungapped Lambda and the database's scale: the product
/// is a double, converted to `Int4` once; `None` where that conversion is undefined (`int4_of`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:216-221 (cited at
/// `RpsCutoffs`)
fn ungapped_x_dropoff(bits: f64, lambda: f64, scale: f64) -> Option<i32> {
    int4_of(scale * (bits * std::f64::consts::LN_2 / lambda).ceil())
}

/// NCBI's `(Int4)` conversion of a double where it is defined: the value truncated toward zero,
/// `None` where that leaves `Int4` (or the value is not a number), which is undefined behaviour in
/// C (C11 6.3.1.4).
fn int4_of(value: f64) -> Option<i32> {
    let truncated = value.trunc();
    (truncated >= f64::from(i32::MIN) && truncated <= f64::from(i32::MAX))
        .then_some(truncated as i32)
}

/// The per-query cutoffs of the preliminary stage at the database's scale.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:216-221
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:340-373
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:937-982
struct RpsCutoffs {
    x_dropoff: i32,
    word_cutoff: i32,
    hit_cutoff: i32,
}

impl RpsRun<'_> {
    /// The query's cutoffs; the refusal where the ungapped X-drop (`-xdrop_ungap`) or a scaled
    /// cutoff leaves `i32`.
    fn cutoffs(
        &self,
        query: &RpsFrame<'_>,
        prelim_space: i64,
        frames: usize,
        scale: f64,
    ) -> Result<RpsCutoffs> {
        let int_scale = scale as i32;
        let x_dropoff = ungapped_x_dropoff(self.x_drop_ungapped_bits, query.kbp_std.lambda, scale)
            .ok_or_else(|| {
                x_drop_overflow(
                    "-xdrop_ungap",
                    self.options.x_drop_ungapped_bits,
                    self.x_drop_ungapped_bits,
                    scale,
                    self.program,
                )
            })?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:926-931
        // (`frames` is 1 for rpsblast, `NUM_FRAMES` for rpstblastn; `prelim_space` the
        // preliminary stage's space, `RpsRun::prelim_spaces`.)
        let hit_cutoff = EvalueStats::new(self.gumbel, prelim_space / frames as i64)
            .hit_cutoff(
                self.prelim_evalue,
                &self.gapped,
                (query.aa_seq.len() - 2) as i32,
                self.cutoff_subject_length,
            )
            .checked_mul(int_scale)
            .ok_or_else(|| scale_overflow(scale, self.program))?;
        let gap_trigger = crate::core::x_drop::gap_trigger_raw_score(
            self.gap_trigger_bits,
            query.kbp_std.lambda,
            query.kbp_std.k,
        );
        Ok(RpsCutoffs {
            x_dropoff,
            word_cutoff: gap_trigger
                .checked_mul(int_scale)
                .ok_or_else(|| scale_overflow(scale, self.program))?
                .min(hit_cutoff),
            hit_cutoff,
        })
    }

    /// The gapped X-drops of the preliminary stage and of the final traceback at the database's
    /// scale (`ext_params->gap_x_dropoff`, `gap_x_dropoff_final`); the refusal where one of
    /// them, or the X-drop of the composition-based redo when it runs
    /// (`gapping_params->x_dropoff`), leaves `i32` (`-xdrop_gap`, `-xdrop_gap_final`). A gapped
    /// X-drop of 0 or below is NCBI's too: the alignment kernels raise it to the gap costs
    /// (`pssm_search::gapalign`, NCBI's `ALIGN_EX`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:455-469
    /// (`min_lambda` is the Lambda of `kbp_gap_std`, BLOSUM62 11/1's, for every query; the
    /// products are `Int4` multiplications.) The redo's X-drop is `RedoRun::x_dropoff_of` of
    /// Lambda over the scaling factor (`s_RescaleSearch`) and the scaled final X-drop
    /// (blast_kappa.c:2375-2383, cited at `RedoRun::new`).
    ///
    /// The kernels then size their memory from the X-drop over the gap extension, in `Int4`
    /// (`kernel_x_drop_fits`); `gap_open` and `gap_extend` are the scaled gap costs.
    fn gapped_x_drops(&self, scale: f64, gap_open: i32, gap_extend: i32) -> Result<(i32, i32)> {
        let int_scale = scale as i32;
        let lambda = self.gapped.lambda;
        let scaled = |raw: i32| {
            if scale > 1.0 {
                raw.checked_mul(int_scale)
            } else {
                Some(raw)
            }
        };
        // An X-drop that the user did not give leaves 32 bits only at an extreme scale: the
        // refusal then names the database's scale.
        let gapped_overflow = || {
            x_drop_overflow(
                "-xdrop_gap",
                self.options.x_drop_gapped_bits,
                self.x_drop_gapped_bits,
                scale,
                self.program,
            )
        };
        let final_overflow = || {
            x_drop_overflow(
                "-xdrop_gap_final",
                self.options.x_drop_final_bits,
                self.x_drop_final_bits,
                scale,
                self.program,
            )
        };
        let raw_gapped = int4_of(self.x_drop_gapped_bits * std::f64::consts::LN_2 / lambda)
            .ok_or_else(gapped_overflow)?;
        let raw_final = int4_of(
            (self.x_drop_final_bits * std::f64::consts::LN_2 / lambda).max(f64::from(raw_gapped)),
        )
        .ok_or_else(final_overflow)?;
        let x_drop_gapped = scaled(raw_gapped).ok_or_else(gapped_overflow)?;
        let x_drop_final = scaled(raw_final).ok_or_else(final_overflow)?;
        // The preliminary stage's kernel runs at every mode; the CBS 0 traceback's (the final
        // X-drop) only without the composition-based redo, whose own gap extension (the costs
        // times 32) keeps its quotient below 2^26.
        let redo = self.sbp.is_some() && self.cbs > 0;
        if !kernel_x_drop_fits(x_drop_gapped, gap_open, gap_extend) {
            return Err(gapped_overflow());
        }
        if !redo && !kernel_x_drop_fits(x_drop_final, gap_open, gap_extend) {
            return Err(final_overflow());
        }
        // The redo converts its X-drop at its own Lambda (`RedoRun::new`);
        // it is checked here, before any search.
        if redo {
            int4_of(RedoRun::x_dropoff_of(
                self.x_drop_final_bits,
                lambda / SCALING_FACTOR,
                x_drop_final,
            ))
            .ok_or_else(final_overflow)?;
        }
        Ok((x_drop_gapped, x_drop_final))
    }

    /// Every query's traced results on the volume `files`, whose profiles' `.psq` letters are
    /// `letters`, in query order, on the scratch of the run's
    /// workers (`workers`): the volume's refusals before the queries' cutoffs (`Err`), else the
    /// queries' results (`RpsVolumeResults`), which a caller reads batch by batch in the order
    /// NCBI's search of a batch meets its errors (the queries search at the same time on several
    /// threads).
    pub(crate) fn search_volume(
        &self,
        files: &RpsVolumeFiles,
        letters: &super::volumes::ProfileLetters<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
        workers: &RpsWorkers,
    ) -> Result<RpsVolumeResults> {
        let rows = files.rows();
        let offsets = files.offsets();
        let scale = files.aux.scale;
        // (A whole number of 1 or more, `rps_hooks::rps_search`.)
        if scale > f64::from(i32::MAX) {
            return Err(scale_overflow(scale, self.program));
        }
        let int_scale = scale as i32;
        let table = RpsTable {
            cells: files.cells().expect("checked when the volume was opened"),
            overflow: files
                .overflow()
                .expect("checked when the volume was opened"),
            pv: files.pv(),
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:558-563
        let (Some(gap_open), Some(gap_extend)) = (
            files.aux.gap_open.checked_mul(int_scale),
            files.aux.gap_extend.checked_mul(int_scale),
        ) else {
            return Err(scale_overflow(scale, self.program));
        };
        // The X-drops at the scale (`gapped_x_drops`, NCBI's `BlastExtensionParametersNew`).
        let (x_drop_gapped, x_drop_final) = self.gapped_x_drops(scale, gap_open, gap_extend)?;
        // Every query's cutoffs, before the search (`BlastHitSavingParametersNew` and
        // `BlastInitialWordParametersNew` set them for the whole batch before its scan); a query
        // whose cutoffs are refused keeps that refusal apart (`RpsVolumeResults::batches` reads it
        // before any search error of its batch) and is not searched.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1372-1374
        let mut refused: Vec<Option<anyhow::Error>> = Vec::with_capacity(self.queries.len());
        let mut cutoffs: Vec<Vec<RpsCutoffs>> = Vec::with_capacity(self.queries.len());
        for (q_idx, query) in self.queries.iter().enumerate() {
            let frames = query.frames.len();
            let frame_cutoffs = |frame: &RpsFrame<'_>, prelim_space: i64| {
                if frame.valid {
                    self.cutoffs(frame, prelim_space, frames, scale)
                } else {
                    Ok(RpsCutoffs {
                        x_dropoff: 0,
                        word_cutoff: 0,
                        hit_cutoff: 0,
                    })
                }
            };
            // The preliminary stage's spaces are the contexts' own (`search_space`) unless the run
            // has its own (`prelim_spaces`).
            let frames = match self.prelim_spaces {
                None => query
                    .frames
                    .iter()
                    .map(|frame| frame_cutoffs(frame, frame.search_space))
                    .collect::<Result<Vec<_>>>(),
                Some(spaces) => query
                    .frames
                    .iter()
                    .zip(&spaces[q_idx * frames..])
                    .map(|(frame, &space)| frame_cutoffs(frame, space))
                    .collect::<Result<Vec<_>>>(),
            };
            match frames {
                Ok(frames) => {
                    cutoffs.push(frames);
                    refused.push(None);
                }
                Err(refusal) => {
                    cutoffs.push(Vec::new());
                    refused.push(Some(refusal));
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1862-1895
        // NCBI divides Lambda by the scale for each HSP's e-value and multiplies it back; at a
        // scale where that does not give Lambda back, NCBI's Lambda changes from HSP to HSP. A
        // refusal of the volume after the queries' cutoffs (NCBI's e-values come after its set-up),
        // so a batch's refusal of its cutoffs comes first; nothing is searched.
        let lambda = self.gapped.lambda;
        if lambda / scale * scale != lambda {
            let searched = refused.iter().map(|_| None).collect();
            return Ok(RpsVolumeResults {
                refused,
                late: Some(anyhow::anyhow!(
                    "a profile database of the scale {scale} (NCBI BLAST+ divides Lambda by the scale for each e-value and multiplies it back, which changes Lambda at this scale) is not supported by LOSAT's {}",
                    self.program
                )),
                searched,
            });
        }
        // Each searched query's result (`None` for a query refused at its cutoffs).
        let mut results: Vec<std::sync::Mutex<Option<Result<RpsTraced>>>> = self
            .queries
            .iter()
            .map(|_| std::sync::Mutex::new(None))
            .collect();
        let redo = self.sbp.filter(|_| self.cbs > 0).map(|sbp| RpsRedoParams {
            files,
            sbp,
            gap_open: files.aux.gap_open,
            gap_extend: files.aux.gap_extend,
            x_drop_final_bits: self.x_drop_final_bits,
            smith_waterman: self.options.use_sw_tback,
            program: self.program,
        });
        let traceback = RpsTraceback {
            rows,
            offsets,
            // Only the traceback without the redo reads it (`RpsTraceback::row_after`).
            row_after_rows: redo
                .is_none()
                .then(|| super::profile::RowAfter::new(files.row_after_rows())),
            letters,
            alphabet: files.alphabet,
            scale,
            karlin_k: &files.aux.karlin_k,
            gap_open,
            gap_extend,
            x_drop_final,
            gapped: self.gapped,
            expect_value: self.evalue,
            hitlist_size: self.hitlist_size,
            redo,
            program: self.program,
        };
        let prelim_karlin = KarlinParams {
            lambda: self.gapped.lambda / scale,
            ..self.gapped
        };
        // `ADAPTIVE_CBS`: the close of the volume's stream trims each query's preliminary hit
        // list (`BlastHSPCBSStreamClose`, blast_traceback.c:1839-1846; `core::blast_hspstream`).
        let cbs_close =
            crate::core::blast_hspstream::adaptive_cbs_close(self.cbs > 0, self.hitlist_size);
        // The concatenated length `Lc` of the volume's unit, NCBI's `query->length` of the
        // search: the length of the diagonal table (`RpsWordScratch::prepare_for`), the bound of
        // the extensions to the right (`RpsRows`, within 0 and the rows) and the end of the
        // interval tree (`RpsPrelim::concat_length`). Where `Lc` is above the rows the extensions
        // read NCBI's row after the rows (`RpsRows::past`, `RpsVolumeFiles::row_after_rows`), and
        // one that reads what NCBI cannot is LOSAT's refusal (`RpsPastRows::overran`).
        let concat_length = self.concat_length;
        let extension_length = concat_length.clamp(0, rows.len() as i32) as usize;
        let past_row =
            (i64::from(concat_length) > rows.len() as i64).then(|| files.row_after_rows());
        // One for each context of each query (query after query, the same number of contexts
        // each), with its flag.
        let past_rows: Option<Vec<RpsPastRows<'_>>> = past_row.as_ref().map(|row| {
            self.queries
                .iter()
                .flat_map(|query| &query.frames)
                .map(|_| RpsPastRows {
                    row,
                    concat_length,
                    overran: std::sync::atomic::AtomicBool::new(false),
                })
                .collect()
        });
        // The preliminary stage of one context of a query (NCBI's frame loop: each context
        // scanned as a subject with its own cutoffs).
        let prelim_of = |q_idx: usize, context: usize| -> RpsPrelim<'_> {
            let frame = &self.queries[q_idx].frames[context];
            let cutoffs = &cutoffs[q_idx][context];
            RpsPrelim {
                word: RpsWordParams {
                    rows: RpsRows {
                        rows,
                        length: extension_length,
                        past: past_rows
                            .as_deref()
                            .map(|past| &past[q_idx * self.queries[q_idx].frames.len() + context]),
                    },
                    table,
                    window: self.window,
                    cutoff: cutoffs.word_cutoff,
                    dropoff: cutoffs.x_dropoff,
                },
                rows,
                offsets,
                concat_length,
                gap_open,
                gap_extend,
                x_drop_gapped,
                hit_cutoff: cutoffs.hit_cutoff,
                karlin: prelim_karlin,
                stats: frame.stats(self.gumbel),
                prelim_evalue: self.prelim_evalue,
                prelim_hitlist_size: self.prelim_hitlist_size,
                query_index: q_idx as u32,
                program: self.program,
            }
        };
        // The units of the volume's search (D2 of PS12s): each valid context of each query whose
        // cutoffs exist, query by query in the context order, each searched on its own
        // (`RpsPrelim::search_frame`) with its HSPs kept apart. A context's search reads nothing
        // another context writes: the diagonal table's entries of earlier contexts read as empty
        // (`wordfinder.rs`) and each context's scan starts its own runs of hits.
        let mut units: Vec<(usize, usize)> = Vec::new();
        // The longest context searched, which bounds the diagonals of the workers' tables
        // (`wordfinder::diag_table_length`).
        let mut longest = 0usize;
        let mut first_unit = Vec::with_capacity(self.queries.len());
        let mut remaining = Vec::with_capacity(self.queries.len());
        for (q_idx, query) in self.queries.iter().enumerate() {
            first_unit.push(units.len());
            let searched = refused[q_idx].is_none();
            if searched {
                for (context, frame) in query.frames.iter().enumerate() {
                    if frame.valid {
                        units.push((q_idx, context));
                        longest = longest.max(frame.length());
                    }
                }
            }
            let count = units.len() - first_unit[q_idx];
            if searched && count == 0 {
                *results[q_idx]
                    .get_mut()
                    .unwrap_or_else(|poisoned| poisoned.into_inner()) = Some(Ok(RpsTraced {
                    hit_list: None,
                    k_mark: None,
                }));
            }
            remaining.push(std::sync::atomic::AtomicUsize::new(count));
        }
        let searched_units: Vec<std::sync::Mutex<Option<Result<Vec<BlastpPreliminaryHsp>>>>> =
            units.iter().map(|_| std::sync::Mutex::new(None)).collect();
        // The rest of a query's search on this volume, run by the worker that searched its last
        // context: its contexts' HSPs concatenated in the context order (`Blast_HSPListAppend`
        // of NCBI's frame loop), the e-values, the reap and the collector (`finish`), then the
        // traceback. A context whose search failed gives the query its error, the first in the
        // context order (where NCBI's loop, and the search before PS12s, stopped).
        let finish_query = |scratch: &mut RpsScratch, q_idx: usize| -> Result<RpsTraced> {
            let query = &self.queries[q_idx];
            let none = RpsTraced {
                hit_list: None,
                k_mark: None,
            };
            let mut combined = std::mem::take(&mut scratch.prelim.combined);
            combined.clear();
            let mut failed = None;
            let mut last = 0;
            for unit in first_unit[q_idx]..first_unit[q_idx] + units_of(&first_unit, &units, q_idx)
            {
                let searched = searched_units[unit]
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner())
                    .take()
                    .expect("each context of a query is searched before its last one finishes");
                match searched {
                    Ok(hsps) if failed.is_none() => combined.extend(hsps),
                    Ok(_) => {}
                    Err(error) => {
                        failed.get_or_insert(error);
                    }
                }
                last = units[unit].1;
            }
            let finished = match failed {
                Some(error) => Err(error),
                None => Ok(prelim_of(q_idx, last).finish(&mut combined, query.stat_length())),
            };
            combined.clear();
            scratch.prelim.combined = combined;
            let Some(mut lists) = finished? else {
                return Ok(none);
            };
            if cbs_close {
                lists.cbs_close_trim(self.hitlist_size);
            }
            let mut traced = traceback.trace_query(
                lists,
                query,
                q_idx as u32,
                self.gumbel,
                &mut scratch.trace,
            )?;
            if let Some(hit_list) = traced.hit_list.as_mut() {
                if query
                    .frames
                    .iter()
                    .any(|frame| frame.aa_seq_nomask.is_some())
                {
                    formatter_identities(hit_list, query, |oid| letters.of(oid))?;
                }
                finish_volume_list(
                    hit_list,
                    self.options.max_hsps_per_subject,
                    self.options.query_cov_hsp_perc,
                    query,
                    self.hitlist_size,
                    self.sorts_results,
                );
            }
            Ok(traced)
        };
        // One unit: its context's HSPs, kept until the query's last context is searched; the
        // worker that searches that one finishes the query at once, so a query's preliminary
        // HSPs live only while the query is in flight.
        let run_unit = |scratch: &mut RpsScratch, unit: usize| {
            let (q_idx, context) = units[unit];
            let frame = &self.queries[q_idx].frames[context];
            // The worker's table (of the unit's length) and buckets (of the volume's rows) for
            // this volume (what the worker kept stays).
            scratch
                .prelim
                .word
                .prepare_for(rows.len(), concat_length, self.window, longest);
            let mut hsps = Vec::new();
            let searched = prelim_of(q_idx, context)
                .search_frame(
                    &mut scratch.prelim,
                    &frame.aa_seq[1..],
                    frame.frame,
                    &mut hsps,
                )
                .map(|()| hsps);
            *searched_units[unit]
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner()) = Some(searched);
            if remaining[q_idx].fetch_sub(1, std::sync::atomic::Ordering::AcqRel) == 1 {
                let traced = out_of_line(|| finish_query(scratch, q_idx));
                *results[q_idx]
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner()) = Some(traced);
            }
        };
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        {
            // For DELTA-BLAST's domain search, bound the number of concurrent diagonal tables so
            // peak RSS on large single-volume domain databases (e.g. `cdd_delta`, 32 MiB per
            // table) stays at or below NCBI's single-threaded `x_FindDomainHits` at 16+ threads.
            const DELTA_DOMAIN_DIAG_BUDGET: usize = 192 * 1024 * 1024;
            let diag_bytes =
                super::wordfinder::diag_table_length(concat_length, self.window, rows.len(), longest)
                    as usize
                    * std::mem::size_of::<u32>();
            let budget_workers = if self.options.domain_search && diag_bytes > 0 {
                (DELTA_DOMAIN_DIAG_BUDGET / diag_bytes).max(1)
            } else {
                usize::MAX
            };
            let active = pool.threads().min(units.len()).min(budget_workers);
            if pool.enabled() && active > 1 {
                use rayon::prelude::*;
                let next = std::sync::atomic::AtomicUsize::new(0);
                pool.install(|| {
                    (0..active).into_par_iter().for_each(|_| {
                        workers.with(|scratch| loop {
                            let unit = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            if unit >= units.len() {
                                break;
                            }
                            run_unit(scratch, unit);
                        });
                    });
                });
                return Ok(take_results(refused, results));
            }
        }
        let _ = pool;
        // Without a pool (or with a single unit) the units run on the caller's thread on one
        // scratch from the free list.
        if !units.is_empty() {
            workers.with(|scratch| {
                for unit in 0..units.len() {
                    run_unit(scratch, unit);
                }
            });
        }
        Ok(take_results(refused, results))
    }
}

// No NCBI counterpart: a call boundary; it computes no value NCBI computes.
/// Calls `run` (the rest of a query's search after its last unit, `RpsRun::search_volume`) in a
/// function of its own that is never inlined, so that the unit's search, whose scan is inlined in
/// the caller, does not take its code: the scan's registers and stack slots then do not depend on
/// the traceback's code. A change there had moved them and made the scan's branches mispredict
/// 3-4 % more (CDD searches 0.9-1.3 % slower; series PR, item PS10e).
#[inline(never)]
fn out_of_line<T>(run: impl FnOnce() -> T) -> T {
    run()
}

/// The number of units of the query `q_idx` (`search_volume`).
fn units_of(first_unit: &[usize], units: &[(usize, usize)], q_idx: usize) -> usize {
    first_unit.get(q_idx + 1).copied().unwrap_or(units.len()) - first_unit[q_idx]
}

/// A volume's search of a run's queries (`RpsRun::search_volume`), read query batch by query
/// batch (`batches`) in the order NCBI's search of one batch meets its errors.
pub(crate) struct RpsVolumeResults {
    /// Each query's refusal of its cutoffs (`None` where they exist).
    refused: Vec<Option<anyhow::Error>>,
    /// The volume's refusal after the queries' cutoffs (the Lambda at the volume's scale): no
    /// query was searched.
    late: Option<anyhow::Error>,
    /// Each searched query's result: its traced lists, or the first error of its search (`None`
    /// for a query that was not searched).
    searched: Vec<Option<Result<RpsTraced>>>,
}

impl RpsVolumeResults {
    /// The results of the query batches whose numbers of queries are `spans`, in query order:
    /// each batch's traced lists, up to the first batch that fails, whose error ends the list.
    /// A batch's error is the first refusal of its queries' cutoffs (NCBI sets the cutoffs of the
    /// whole batch before its scan, `search_volume`), else the volume's refusal after them, else
    /// the first error of its queries' searches in query order (the queries search at the same
    /// time).
    pub(crate) fn batches(self, spans: &[usize]) -> Vec<Result<Vec<RpsTraced>>> {
        let Self {
            refused,
            mut late,
            searched,
        } = self;
        let mut refused = refused.into_iter();
        let mut searched = searched.into_iter();
        let mut batches = Vec::with_capacity(spans.len());
        for &count in spans {
            // The batch's refusals are all consumed: the next batch never reads them.
            let first_refusal = refused
                .by_ref()
                .take(count)
                .fold(None, |first, refusal| first.or(refusal));
            let results: Vec<Option<Result<RpsTraced>>> = searched.by_ref().take(count).collect();
            let batch = match first_refusal.or_else(|| late.take()) {
                Some(refusal) => Err(refusal),
                None => results
                    .into_iter()
                    .map(|result| result.expect("a query whose cutoffs exist is searched"))
                    .collect(),
            };
            let failed = batch.is_err();
            batches.push(batch);
            if failed {
                break;
            }
        }
        batches
    }

    /// The results of the run's queries as one batch (`batches`), for the tests of the order of
    /// the errors.
    #[cfg(test)]
    pub(crate) fn all(self) -> Result<Vec<RpsTraced>> {
        let count = self.refused.len();
        self.batches(&[count])
            .pop()
            .expect("one batch has its result")
    }
}

/// Every query's result of a volume's search, in query order (`RpsVolumeResults`).
fn take_results(
    refused: Vec<Option<anyhow::Error>>,
    results: Vec<std::sync::Mutex<Option<Result<RpsTraced>>>>,
) -> RpsVolumeResults {
    let searched = results
        .into_iter()
        .zip(&refused)
        .map(|(result, refusal)| {
            let result = result
                .into_inner()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            match refusal {
                Some(_) => None,
                None => Some(result.expect("every searched query of the volume has its result")),
            }
        })
        .collect();
    RpsVolumeResults {
        refused,
        late: None,
        searched,
    }
}

/// The scratch of one worker of an RPS search: its preliminary stage's (the word finder's
/// diagonal table and buckets, the gapped kernels' memory) and its traceback's.
struct RpsScratch {
    prelim: RpsPrelimScratch,
    trace: RpsTraceScratch,
}

impl RpsScratch {
    /// A worker's scratch, its buckets' blocks from the run's `pool`.
    fn new(pool: &std::sync::Arc<RpsBlockPool>) -> Self {
        Self {
            prelim: RpsPrelimScratch::new(RpsWordScratch::new(pool.clone())),
            trace: RpsTraceScratch::new(),
        }
    }
}

/// The scratch of a run's workers (series PR, item PS12s, D3): handed out from a free list so
/// the number of live diagonal tables never exceeds the peak number of concurrently running units
/// (e.g. 1 for a single-query batch at any thread count), kept from volume to volume inside one
/// search where each volume would otherwise allocate and fault in its own tables. A worker's
/// diagonal table stays where the next volume's table has the same length
/// (`RpsWordScratch::prepare_for`); the buckets' blocks come from one pool of the run
/// (`RpsBlockPool`, M2), to which each bucket gives its blocks back as soon as its hits are read.
/// A search's values do not depend on what its scratch held before (`wordfinder.rs`), so only the
/// memory and its reuse change.
pub(crate) struct RpsWorkers {
    free: std::sync::Mutex<Vec<RpsScratch>>,
    pool: std::sync::Arc<RpsBlockPool>,
}

impl RpsWorkers {
    /// An empty free list of worker scratches (allocated lazily up to peak concurrent units).
    pub(crate) fn new(_threads: usize) -> Self {
        Self {
            free: std::sync::Mutex::new(Vec::new()),
            pool: std::sync::Arc::new(RpsBlockPool::default()),
        }
    }

    fn take(&self) -> RpsScratch {
        self.free
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .pop()
            .unwrap_or_else(|| RpsScratch::new(&self.pool))
    }

    fn put(&self, scratch: RpsScratch) {
        self.free
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .push(scratch);
    }

    /// Runs `work` on a scratch taken from the free list, returning it when `work` finishes.
    fn with<T>(&self, work: impl FnOnce(&mut RpsScratch) -> T) -> T {
        struct ReturnGuard<'a> {
            workers: &'a RpsWorkers,
            scratch: Option<RpsScratch>,
        }
        impl Drop for ReturnGuard<'_> {
            fn drop(&mut self) {
                if let Some(scratch) = self.scratch.take() {
                    self.workers.put(scratch);
                }
            }
        }
        let mut guard = ReturnGuard {
            workers: self,
            scratch: Some(self.take()),
        };
        work(guard.scratch.as_mut().expect("held until drop"))
    }
}

/// What the RPS searches of one database keep from search to search (series PR, item PS13a:
/// deltablast's domain search searches each query batch of its run, `blast_engine::delta_hooks`):
/// the set-up of the volumes, which checks each volume's files at its first search only
/// (`RpsVolumeSetup::open`). Worker scratches (`RpsWorkers`) live only for the duration of each
/// RPS search and are dropped before the main `-db` search or report prefetch, matching NCBI's
/// `x_FindDomainHits` lifetime. The volumes' mappings are each search's own, as NCBI maps the
/// files for each search (`CSetupFactory::CreateRpsStructures`): their pages leave the resident
/// set when the search ends.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:72-83
#[derive(Default)]
pub(crate) struct RpsRunState {
    pub(crate) setup: Option<super::volumes::RpsVolumeSetup>,
}

/// The gapped K that the outfmt 0 footer of each query prints: for the first query of a batch,
/// `RPS_K_MULT` times the `.aux` K of the last list the batch's traceback reads, else the
/// table's K; every other query prints the table's K. Where the stream sorts by score
/// (`by_score`: rpsblast with composition-based statistics, `stream_lists`) the lists are read
/// query by query in the order of `s_EvalueCompareHSPLists` and each writes the K (the mark is
/// that of the last query that has one); else the lists of the batch are read by ordinal id
/// (the mark is the highest): without composition-based statistics only those whose rescaled
/// PSSM exists write it, and with them (rpstblastn's redo) each does.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1234-1235
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:292-298
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:98-103
/// The traceback writes only `kbp_gap[0]`, the first query's (`PR/STATE.md` point 1, checked on
/// the survey's 39 outfmt 0 reports); an invalid query has no Karlin block in its footer.
pub(crate) fn batch_gapped_k(
    batches: &[std::ops::Range<usize>],
    k_marks: &[Option<u32>],
    karlin_k: &[f64],
    table_k: f64,
    by_score: bool,
) -> Vec<f64> {
    let mut gapped_k = vec![table_k; k_marks.len()];
    for batch in batches {
        if batch.is_empty() {
            continue;
        }
        let marks = &k_marks[batch.clone()];
        let mark = if by_score {
            marks.iter().flatten().next_back()
        } else {
            marks.iter().flatten().max()
        };
        if let Some(mark) = mark {
            gapped_k[batch.start] = RPS_K_MULT * karlin_k[*mark as usize];
        }
    }
    gapped_k
}

/// A volume's final hit list of a query: the filters (`-max_hsps`, `-qcov_hsp_perc`), then the
/// list sorted by E-value where the unit's total length is above 0 (`sorts_results`; the
/// traceback sorted it before the filters, `traceback::RpsTraceback::trace_query`,
/// blast_traceback.c:1306) and cut to the hit list size, before the volumes are merged.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1777
fn finish_volume_list(
    hit_list: &mut crate::algorithm::blastp::hsp::BlastpHitList,
    max_hsps: usize,
    query_cov_hsp_perc: f64,
    query: &RpsQuery<'_>,
    hitlist_size: usize,
    sorts_results: bool,
) {
    filter_results(hit_list, max_hsps, query_cov_hsp_perc, |frame| {
        query.context(frame).first_query_length
    });
    if sorts_results {
        hit_list.sort_by_evalue();
    }
    hit_list.prune_by_size(hitlist_size);
}

#[cfg(test)]
#[path = "engine_tests.rs"]
mod tests;
