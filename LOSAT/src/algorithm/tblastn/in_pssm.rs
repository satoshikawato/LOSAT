//! tblastn `-in_pssm F` (PSI-TBLASTN; line PS, item PS7, DECISIONS 35): the PSSM read as psiblast
//! reads it (`algorithm::psiblast::in_pssm`, with tblastn's program name and matrix), its query
//! record, the lookup table of its rows, the program name of its reports, and LOSAT's refusals of
//! what NCBI BLAST+ 2.17.0's tblastn cannot search with a PSSM.
//!
//! The search is tblastn's (the translation of the subjects, the units, chunks and threads, the
//! containment trees, the fence retry, the `-db` stream, the hit lists and the reports): only the
//! source of the scores (`score::PssmScore`), the Karlin-Altschul blocks and the lookup table are
//! the PSSM's, as NCBI's `eBlastTypePsiTblastn` (`PSSM_QUERY_MASK | eBlastTypeTblastn`) takes
//! tblastn's code everywhere except where `Blast_QueryIsPssm` selects the PSSM.
//!
//! NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_program.h:82-83
//! ```c
//!     /** PSI Tblastn */
//!     eBlastTypePsiTblastn = (PSSM_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK | TRANSLATED_SUBJECT_MASK),
//! ```

use std::path::Path;

use anyhow::{bail, Result};

use super::args::{ResolvedTblastnArgs, TblastnArgs};
use crate::algorithm::psiblast::in_pssm::{read_in_pssm, PssmReader};
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::algorithm::pssm_search::PsiSearch;
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;
use crate::blastinput::fasta_reader::FastaRecord;

/// The program line of PSI-TBLASTN's reports (`Blast_ProgramNameFromType` of
/// `eBlastTypePsiTblastn`, upper-cased by the formatter).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:690
/// ```c++
///     case eBlastTypePsiTblastn:  retval.assign("psitblastn"); break;
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:94
/// ```c++
///     m_Program = Blast_ProgramNameFromType(opts.GetProgramType());
/// ```
pub(super) const PROGRAM: &str = "PSITBLASTN";

/// Reads `-in_pssm path` for tblastn: psiblast's reader (NCBI's `CPsiBlastArgs`, the same handler
/// for both programs), whose refusals name LOSAT's TBLASTN and whose scores of a PSSM of
/// frequency ratios only are made with the options' matrix (`matrix_name`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:128-130
/// ```c++
///     m_PsiBlastArgs.Reset(new CPsiBlastArgs(CPsiBlastArgs::eNucleotideDb));
///     arg.Reset(m_PsiBlastArgs);
///     m_Args.push_back(arg);
/// ```
pub(super) fn read(path: &Path, matrix_name: &str) -> Result<PsiSearch> {
    let (search, _file) = read_in_pssm(
        path,
        PssmReader {
            program: "TBLASTN",
            matrix_name,
        },
    )?;
    Ok(search)
}

/// The check of `-db_gencode` with `-in_pssm`: the value is not used (PSI-TBLASTN translates with
/// code 1, `check_options`), and NCBI's command line takes the codes of its list only, without
/// the code 32 that LOSAT's TBLASTN adds for its own translations (`PD-TLOSAN-LOCAL-GENCODE-32`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1011-1015
/// ```c++
///     static int gcs[] = {1,2,3,4,5,6,9,10,11,12,13,14,15,16,
///                         21,22,23,24,25,26,27,28,29,30,31,33};
/// ```
pub(super) fn check_db_gencode(code: u8) -> Result<()> {
    if code == 32 {
        bail!(
            "-db_gencode 32 with -in_pssm is not supported by LOSAT's TBLASTN (NCBI BLAST+ 2.17.0 does not take the code 32, and its PSI-TBLASTN translates the subjects with code 1)"
        );
    }
    Ok(())
}

/// The refusal of `-in_pssm` where the queries are given as records (`run_local`,
/// `run_local_database`: the library API refuses the option first,
/// `api::search::pssm_file_options`, and the Web's checks, `args::check_options`).
pub(super) fn check_records_only(args: &TblastnArgs) -> Result<()> {
    if args.in_pssm.is_some() {
        bail!("-in_pssm (a PSSM query) is not supported by LOSAT's TBLASTN where the queries are given as records");
    }
    Ok(())
}

/// The refusal of what NCBI BLAST+ 2.17.0's PSI-TBLASTN runs and LOSAT's does not.
fn refuse(what: impl std::fmt::Display) -> anyhow::Error {
    anyhow::anyhow!("{what} with -in_pssm is not supported by LOSAT's TBLASTN")
}

/// LOSAT's refusals of a PSI-TBLASTN search, after NCBI's checks of the options and the opening of
/// the database, before anything is written: NCBI's crash of `-mt_mode 1` (DECISIONS 35 Q4), the
/// output formats other than 0, 6, 7 and 10 (item PS15b of line PS; 11 crashes NCBI), and a PSSM
/// query long enough for NCBI to split it.
pub(super) fn check_limits(
    args: &TblastnArgs,
    resolved: &ResolvedTblastnArgs,
    psi: &PsiSearch,
    outfmt: i32,
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:159-162
    // ```c++
    // 	if ((mt_mode == CMTArgs::eSplitByQueries) && (num_threads > 1)){
    // 		m_UsageReport.AddParam(CBlastUsageReport::eMTMode, CMTArgs::eSplitByQueries);
    // 		return x_RunMTBySplitQuery();
    // 	}
    // ```
    // NCBI's split by queries reads the query file, which a PSSM search has not opened: NCBI
    // BLAST+ 2.17.0 aborts (probe `ts1.mt1.n2.6`). With `-subject` it runs on one thread.
    if args.mt_mode == Some(1) && resolved.num_threads > 1 && args.db.is_some() {
        return Err(refuse(
            "-mt_mode 1 with -num_threads above 1 and -db, with which NCBI BLAST+ 2.17.0's tblastn crashes,",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:791-796 (the tabular
    // formats of tblastn; the pairwise report, `PrintOneResultSet`)
    // ```c++
    //     if (m_FormatType >= CFormattingArgs::eTabular &&
    //         m_FormatType <= CFormattingArgs::eCommaSeparatedValues) {
    //         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, m_CustomDelim);
    // ```
    // The pairwise report and the tabular formats 6, 7 and 10 (their field lists too).
    if !matches!(outfmt, 0 | 6 | 7 | 10) {
        bail!(
            "output format {outfmt} with -in_pssm is not supported by LOSAT's TBLASTN yet{}",
            if outfmt == 11 {
                " (NCBI BLAST+ 2.17.0's tblastn stops with a C++ exception)"
            } else {
                ""
            }
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:91-94
    // ```c++
    //     else if (Blast_ProgramIsPsiBlast(program))
    //         retval = 10000;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:191-196 (the
    // chunks of a split query are set up without the PSSM)
    // A PSSM query that NCBI splits is searched with the matrix in its preliminary stage, as
    // psiblast's (`algorithm::psiblast::in_pssm::check_in_pssm_options`).
    let sizes =
        crate::algorithm::blastp::protein_query_split::protein_query_split_sizes("TBLASTN", 10000)?;
    let length = psi.query.seq().len();
    if crate::algorithm::blastp::protein_query_split::split_protein_batch(&[length], sizes)
        .is_some()
    {
        bail!(
            "a PSSM query of {length} letters, which NCBI BLAST+ splits into chunks searched without the PSSM, is not supported by LOSAT's TBLASTN"
        );
    }
    Ok(())
}

/// The query record of a PSSM search: psiblast's record of the PSSM's Bioseq (its Seq-id string,
/// its title, its letters). With `-parse_deflines` NCBI's reports show the Seq-id before the title
/// (`Query= Query_1 <title>`, `qseqid` `Query_1`; a PSSM without a title shows its id alone), as
/// the title of a record whose first word is that id shows them.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:732-734
/// ```c++
///     string all_id_str = GetSeqIdString(cbs, believe_query);
///     all_id_str += " ";
///     all_id_str = NStr::TruncateSpaces(all_id_str + GetSeqDescrString(cbs));
/// ```
pub(super) fn query_record(psi: &PsiSearch, parse_deflines: bool) -> FastaRecord {
    let query = &psi.query;
    if !parse_deflines {
        return query.clone();
    }
    let mut title = query.local_id.clone().into_bytes();
    if !query.title.is_empty() {
        title.push(b' ');
        title.extend_from_slice(&query.title);
    }
    FastaRecord::new(query.local_id.clone(), &title, query.seq())
}

/// The lookup segments of a query of `length` letters whose filtering locations are `masks`
/// (half-open, sorted, not overlapping): NCBI's complement for one valid protein context at offset
/// 0, as inclusive ranges (`left > right` when a mask reaches the end).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1044-1110
/// ```c
///       start_offset = query_info->contexts[context].query_offset;
///       end_offset = query_info->contexts[context].query_length
///           + start_offset - 1;
///       ...
///       if (mask_loc == NULL || mask_loc->seqloc_array[context] == NULL) {
///          tail = BlastSeqLocNew(tail ? &tail : complement_mask,
///                                start_offset, end_offset);
///          continue;
///       }
///       ...
///          if (first) {
///             last_interval_open = TRUE;
///             first = FALSE;
///
///             if (filter_start > start_offset) {
///                /* beginning of sequence not filtered */
///                left = start_offset;
///             } else {
///                /* beginning of sequence filtered */
///                left = filter_end + 1;
///                continue;
///             }
///          }
///
///          right = filter_start - 1;
///
///          tail = BlastSeqLocNew((tail ? &tail : complement_mask), left, right);
///          if (filter_end >= end_offset) {
///             /* last masked region at end of sequence */
///             last_interval_open = FALSE;
///             break;
///          } else {
///             left = filter_end + 1;
///          }
///       }
///
///       if (last_interval_open) {
///          /* Need to finish SSeqRange* for last interval. */
///          right = end_offset;
///          tail = BlastSeqLocNew((tail ? &tail : complement_mask), left, right);
///       }
/// ```
fn lookup_segments(length: usize, masks: &[(usize, usize)]) -> Vec<(i64, i64)> {
    let end_offset = length as i64 - 1;
    if masks.is_empty() {
        return vec![(0, end_offset)];
    }
    let mut segments = Vec::new();
    let mut first = true;
    let mut last_interval_open = true;
    let mut left = 0i64;
    for &(start, end) in masks {
        let (filter_start, filter_end) = (start as i64, end as i64 - 1);
        if first {
            last_interval_open = true;
            first = false;
            if filter_start > 0 {
                left = 0;
            } else {
                left = filter_end + 1;
                continue;
            }
        }
        segments.push((left, filter_start - 1));
        if filter_end >= end_offset {
            last_interval_open = false;
            break;
        }
        left = filter_end + 1;
    }
    if last_interval_open {
        segments.push((left, end_offset));
    }
    segments
}

/// The lookup table of a PSSM query: the neighbouring words of the PSSM's rows over the lookup
/// segments, the complement of the query's filtering locations (`-seg`'s masks, also with
/// `-soft_masking true`, which keeps them out of the table only). `query` is the query record's
/// letters.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:631-637
/// ```c
///     if (lookup_segments) {
///         if (program_number == eBlastTypeBlastx || program_number == eBlastTypeTblastx) {
///         ...
///         } else {
///             BLAST_ComplementMaskLocations(program_number, query_info,
///                                           filter_maskloc, lookup_segments);
///         }
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:86-99 (quoted at
/// `pssm_search::lookup`: a score block with a PSSM indexes the PSSM's rows)
pub(super) fn lookup_table(
    rows: &QueryPssm,
    query: &[u8],
    seg: Option<&crate::utils::seg::SegParams>,
    mask_lowercase: bool,
    word_length: usize,
    threshold: i32,
) -> BlastAaLookupTable {
    let masks =
        super::search_seed::encode_tblastn_lookup_query(query, seg, mask_lowercase).seg_masks;
    let segments = lookup_segments(query.len(), &masks);
    crate::algorithm::pssm_search::lookup::build_pssm_lookup_in(
        rows,
        word_length,
        threshold,
        &segments,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI blast_filter.c:1044-1110: the complement of the masks of one context.
    #[test]
    fn lookup_segments_are_the_complement_of_the_masks() {
        assert_eq!(lookup_segments(10, &[]), vec![(0, 9)]);
        assert_eq!(lookup_segments(10, &[(2, 4)]), vec![(0, 1), (4, 9)]);
        assert_eq!(lookup_segments(10, &[(0, 3)]), vec![(3, 9)]);
        assert_eq!(lookup_segments(10, &[(6, 10)]), vec![(0, 5)]);
        assert_eq!(
            lookup_segments(12, &[(0, 2), (5, 7), (9, 12)]),
            vec![(2, 4), (7, 8)]
        );
        // A mask of the whole query leaves the inverted interval NCBI makes (no word in it).
        assert_eq!(lookup_segments(5, &[(0, 5)]), vec![(5, 4)]);
    }

    // NCBI blast_args.cpp:1011-1015: the code 32 is not NCBI's.
    #[test]
    fn the_genetic_code_32_is_refused_with_a_pssm() {
        assert!(check_db_gencode(1).is_ok());
        assert!(check_db_gencode(33).is_ok());
        assert!(check_db_gencode(32)
            .unwrap_err()
            .to_string()
            .contains("-db_gencode 32 with -in_pssm"));
    }
}
